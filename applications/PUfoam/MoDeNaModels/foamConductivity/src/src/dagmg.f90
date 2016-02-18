! This file is part of AGMG software package,
! Release 3.1.1 built on "Nov 23 2011"
!
!    AGMG is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AGMG is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AGMG.  If not, see <http://www.gnu.org/licenses/>.
!
! Up-to-date copies of the AGMG package can be obtained
! from the Web pages <http://homepages.ulb.ac.be/~ynotay/AGMG>.
!
! You can acknowledge, citing references [1] [2], and [3], the contribution
! of this package in any scientific publication dependent upon it use.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010
!
! [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
!    convergence rate, Report GANMN 10-03, Universite Libre de Bruxelles,
!    Brussels, Belgium, 2010.
!
! [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
!    equations, Report GANMN 11-01, Universite Libre de Bruxelles, Brussels,
!    Belgium, 2011.
!
! See the accompanying userguide for more details on how to use the software,
! and the README file for installation instructions.
!
! AGMG Copyright (C) 2011 Yvan NOTAY
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file provides dagmg and dependencies; dagmg is a sequential
! implementation for real matrices in double precision
! of the method presented in [1], where the used algorithms are described
! in detail. From realease 3.x, the coarsening has been modified according
! to the results in [2,3], see the release notes in the README file
! for more details.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! SUBROUTINE dagmg (MAIN DRIVER): see bottom of this file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! PARAMETERS DEFINITON -  may be tuned by expert users
!-----------------------------------------------------------------------
  MODULE dagmg_mem
    SAVE
!
! INTEGER
!
!  maxlev  is the maximal number of levels
!          (should be large enough - much larger than needed is armless).
!  real_len is the length of 1 REAL(kind(0.0d0)) in byte
!        (used only to display information on memory usage).
!  nsmooth  is the number of pre- and post-smoothing steps;
!  smoothtype indicates which smoother use:
!    if smoothtype==1, the smoothing scheme is
!        pre-smoothing: Forward Gauss-Seidel, then Backward GS, then Fwd GS,etc
!       post-smoothing: Backward Gauss-Seidel, then Forward GS, then Bwd GS,etc
!                       (nsmooth sweeps for both pre- and post-smoothing)
!  nstep   is the maximum number of coarsening steps
!          nstep<0 means that coarsening is stopped only according to
!          the matrix size, see parameter maxcoarsesize.
!  nlvcyc  is the number of coarse levels (from bottom) on which V-cycle
!          formulation is enforced (Rmk: K-cycle always allowed at first
!          coarse level).
!  npass   is the maximal number of pairwise aggregation passes for each
!          coarsening step, according to the algorithms in [2,3].
!  maxcoarsesize: when the size of the coarse grid matrix is less than or
!                 equal to maxcoarsesize*N^(1/3),  it is factorized exactly
!                 and coarsening is stopped;
!         maxcoarsesizeslow: in case of slow coarsening,
!                 exact factorization can be performed when the size of
!                 the coarse grid matrix is less than or equal to
!                 maxcoarsesizeslow*N^(1/3).
!         (N is the number of rows of the input matrix)
    INTEGER, PARAMETER :: maxlev=40, real_len=8
    INTEGER, PARAMETER :: nsmooth=1, smoothtype=1, nstep=-1, nlvcyc=0
    INTEGER, PARAMETER :: npass=2,maxcoarsesize=400,maxcoarsesizeslow=4000
!
! REAL
!
!  resi is the threshold t for the relative residual error in inner FCG & GCR
!       iterations, see Algorithm 3.2 in [1]
!  trspos is a threshold: if a row has a positive offdiagonal entry larger
!         than trspos times the diagonal entry, the corresponding node is
!         transferred unaggregated to the coarse grid
!  kaptg_ is the threshold used to accept or not a tentative aggregate
!         when applying the coarsening algorithms from [2,3];
!         kaptg_blocdia is used for control based on bloc diagonal smoother [2];
!         kaptg_dampJac is used for control based on Jacobi smoother [3].
!  checkdd is the threshold to keep outside aggregation nodes where
!         the matrix is strongly diagonally dominant (based on mean of row
!         and column);
!         In fact, AGMG uses the maximum of |checkdd| and of the default value
!            according to kaptg_ as indicated in [2,3]
!            (hence |checkdd| < 1 ensures that one uses this default value)
!         checkdd <0 : consider |checkdd|, but base the test on minus the
!                sum of offdiagonal elements, without taking the absolute value.
!  targetcoarsefac is the target coarsening factor (parameter tau in the main
!         coarsening algorithm in [2,3]): further pairwise aggregation passes
!         are omitted once the number of nonzero entries has been reduced by a
!         factor of at least targetcoarsefac.
!  fracnegrcsum: if, at some level, more than fracnegrcsum*nl nodes,
!         where nl is the total number of nodes at that level, have
!         negative mean row and column sum, then the aggregation algorithm
!         of [2,3] is modified, exchanging all diagonal entries for the mean
!         row and column sum (that is, the algorithm is applied to a
!         modified matrix with mean row and colum sum enforced to be zero).
    REAL(kind(0.0d0)), PARAMETER :: resi=0.2, trspos=0.45
    REAL(kind(0.0d0)), PARAMETER :: kaptg_blocdia=8, kaptg_dampJac=10
    REAL(kind(0.0d0)), PARAMETER :: checkdd=0.5
    REAL(kind(0.0d0)), PARAMETER :: targetcoarsefac=2.0**npass
    REAL(kind(0.0d0)), PARAMETER :: fracnegrcsum=0.25
!!!!!!!!!!!!!!!!!!!! END of PARAMETERS DEFINITION -----------------
!!!!!!!!!!!!!!!!!!! Internal variables declaration
!
    TYPE InnerData
       REAL(kind(0.0d0)), DIMENSION(:), POINTER :: a
       INTEGER, DIMENSION(:), POINTER :: ja
       INTEGER, DIMENSION(:), POINTER :: ia
       INTEGER, DIMENSION(:), POINTER :: il
       INTEGER, DIMENSION(:), POINTER :: iu
       REAL(kind(0.0d0)), DIMENSION(:), POINTER :: p
       INTEGER, DIMENSION(:), POINTER :: idiag
       INTEGER, DIMENSION(:), POINTER :: ind
       INTEGER, DIMENSION(:), POINTER :: iext
       INTEGER, DIMENSION(:), POINTER :: ilstout
       INTEGER, DIMENSION(:), POINTER :: lstout
       INTEGER, DIMENSION(:), POINTER :: ilstin
       INTEGER, DIMENSION(:), POINTER :: lstin
    END TYPE InnerData
!
    TYPE(InnerData) :: dt(maxlev)
    REAL(kind(0.0d0)), ALLOCATABLE :: scald(:)
    INTEGER :: nn(maxlev),kstat(2,maxlev)=0,innermax(maxlev)
    INTEGER :: nlev,nwrkcum,iout,nrst
    INTEGER :: maxcoarset=maxcoarsesize,maxcoarseslowt=maxcoarsesizeslow
    REAL(kind(0.0d0)) :: memi=0.0,memax=0.0,memr=0.0,mritr,rlenilen
    REAL(kind(0.0d0)) :: wcplex(4),fracnz(maxlev)
    LOGICAL :: spd,wfo,wff,allzero,trans,transint,zerors,gcrin
    REAL(kind(0.0d0)), PARAMETER :: cplxmax=3.0, xsi=0.6d0
    REAL(kind(0.0d0)), PARAMETER :: repsmach=1.49d-8 !SQRT(EPSILON(1.0d0)) !because of absoft
    INTEGER :: nlc(2),nlcp(2),nlc1(2),icum,imult
    REAL(kind(0.0d0)) :: ngl(2),nglp(2),nlctot(2),ngl1(2),ngltot(2)
    INTEGER, DIMENSION(:), POINTER :: iextL1
    INTEGER, PARAMETER :: IRANK=-9999
    REAL(kind(0.0d0)), PARAMETER ::    &
               checkddJ=1.25d0 !MAX(ABS(checkdd),kaptg_dampJac/(kaptg_dampJac-2)) !because of absoft
    REAL(kind(0.0d0)), PARAMETER ::    &
               checkddB=1.28571428571429d0 !MAX(ABS(checkdd),(kaptg_blocdia+1)/(kaptg_blocdia-1)) !because of absoft
  END MODULE dagmg_mem
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of Internal variables declaration
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! TIMING
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_mestime(id,cputm,eltm)
    IMPLICIT NONE
    INTEGER, SAVE :: cpt_init(10)=-1,cpt_fin,cpt_max,freq,cpt
    REAL, SAVE :: t1(10), t2
    REAL(kind(0.0d0)) :: cputm,eltm
    INTEGER :: id
    IF (id>0) THEN
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t2)
       CALL SYSTEM_CLOCK(cpt_fin,freq,cpt_max)
       !
       cpt = cpt_fin - cpt_init(id)
       IF (cpt_fin < cpt_init(id)) cpt = cpt + cpt_max
       eltm = dble(cpt) / freq
       cputm = dble(t2 - t1(id))
       !
    ELSE
       !
       CALL SYSTEM_CLOCK(cpt_init(-id),freq,cpt_max)
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t1(-id))
       !
    END IF
    RETURN
  END SUBROUTINE dagmg_mestime
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of TIMING
!-----------------------------------------------------------------------
MODULE dagmg_ALLROUTINES
CONTAINS
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! MEMORY RELEASE
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_relmem
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l
    !
    DO l=1,nlev-1
       IF(nn(l) > 0) DEALLOCATE(dt(l)%a,dt(l)%ja,dt(l)%il,dt(l)%iu)
       IF (nn(l+1) > 0) DEALLOCATE(dt(l)%ind)
    END DO
    memi=0
    memr=0
    memax=0
    RETURN
  END SUBROUTINE dagmg_relmem
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of MEMORY RELEASE
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! APPLY MG PRECONDITIONER
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_applyprec( N,f,X,a,ja,ia,flop)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER       :: N
    REAL(kind(0.0d0)) ::  flop
    REAL(kind(0.0d0)) ::  f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)), ALLOCATABLE :: S(:)
    !
    mritr=nwrkcum+N
    ALLOCATE (S(nwrkcum+N))
    CALL dagmg_prec_matv(1,f,X,S,flop,S(N+1),.FALSE.)
    kstat(2,1)=kstat(2,1)+1
    DEALLOCATE (S)
    !
    RETURN
  END SUBROUTINE dagmg_applyprec
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of APPLY MG PRECONDITIONER
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! PARTIAL ORDERING OF MATRIX ROWS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_partroword(n, a, ja, ia, idiag, w, iw, iext)
!
!  Performs partial ordering of matrix rows:
!   first lower part, next diagonal, eventually upper part;
!   idiag(i) is set such that it points to diagonal element
!
!  w,iw: working arrays; length should be at least equal to the maximum
!                        number of entries in (the uper part of) a row
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER, TARGET :: iw(*)
    REAL(kind(0.0d0)) :: a(*), p
    REAL(kind(0.0d0)), TARGET :: w(*)
    INTEGER :: i, j, k, ipos, nzu, nze
    DO i = 1, n
       !       find diag, save upper part, concatenate lower part
       ipos = ia(i)
       nzu = 0
       nze = 0
       DO k = ia(i), ia(i+1) - 1
          j = ja(k)
          IF (j.GT.i) THEN
             nzu = nzu + 1
             w(nzu) = a(k)
             iw(nzu) = j
          ELSEIF (j.LT.i) THEN
             a(ipos) = a(k)
             ja(ipos) = j
             ipos = ipos + 1
          ELSE
             p = a(k)
          ENDIF
       ENDDO
       !
       !       copy back diagonal entry
       idiag(i) = ipos
       a(ipos) = p
       ja(ipos) = i
       !
       !       copy back upper part
       DO k = 1, nzu
          ipos = ipos + 1
          a(ipos) = w(k)
          ja(ipos) = iw(k)
       ENDDO
       !
    ENDDO
    RETURN
  END SUBROUTINE dagmg_partroword
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of PARTIAL ORDERING OF MATRIX ROWS
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! CSR TO DLU FORMAT
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrdlu(n,a,ja,ia,idiag,ao,jao,il,iu,trans,iext,iexto)
!
!  Convert matrix in csr format with partial orederings of the rows
!   into separate d+l+u format:
!     first diagonal elements in ao(1:n),
!     next rows of lower part (pointer: il),
!     and thereafter the rows of upper part (pointer: iu)
!
!   If trans==.TRUE. , performs simultaneously the transposition of the matrix
!
!   If trans==.FASLE.,
!     il, iu (iexto) may point to the same memory location as ia, idiag (iext)
!     (if idiag (iext) has length at least n+1 in the calling program)
!
!   Pay attention that joa is defined as jao(n+1:*)
!     (entries jao(1:n) not needed with this format)
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n), jao(n+1:*), il(n+1), iu(n+1)
    INTEGER, OPTIONAL :: iext(n),iexto(n+1)
    REAL(kind(0.0d0)) :: a(*),ao(*)
    INTEGER :: i,j,ili,iui,iei,ipos,ili0,iui0,iei0,nl,nu,k,next
    LOGICAL :: trans
    !
    ! First set il(i), iu(i) and iext(i) equal to the number of elements in
    !       repsective rows;
    !       compute total number of elements in lower and upper part
    !
    nl=0
    nu=0
    IF (.NOT.trans) THEN
       DO i=1,n
          ili=idiag(i)-ia(i)
          iui=ia(i+1)-idiag(i)-1
          iu(i)=iui
          il(i)=ili
          nl=nl+ili
          nu=nu+iui
       END DO
    ELSE
       il(2:n+1)=0
       iu(2:n+1)=0
       DO i=1,n
          iui=idiag(i)-ia(i)
          ili=ia(i+1)-idiag(i)-1
          nl=nl+ili
          nu=nu+iui
          DO k=ia(i),idiag(i)-1
             iu(ja(k)+1)=iu(ja(k)+1)+1
          END DO
          DO k=idiag(i)+1,ia(i+1)-1
             il(ja(k)+1)=il(ja(k)+1)+1
          END DO
       END DO
    END IF
    !
    ! Now perform actual copying in new format
    !
    ipos=0
    ili=n+1
    iui=ili+nl
    iei=iui+nu
    IF (.NOT.trans) THEN
       DO i=1,n
          ili0=ili
          iui0=iui
          DO j=1,il(i)
             ipos=ipos+1
             ao(ili)=a(ipos)
             jao(ili)=ja(ipos)
             ili=ili+1
          END DO
          ipos=ipos+1
          ao(i)=a(ipos)
          DO j=1,iu(i)
             ipos=ipos+1
             ao(iui)=a(ipos)
             jao(iui)=ja(ipos)
             iui=iui+1
          END DO
          iu(i)=iui0
          il(i)=ili0
       END DO
       iu(n+1)=iui
       il(n+1)=ili
    ELSE
       ! compute pointers from lengths
       il(1)=ili
       iu(1)=iui
       DO i=1,n
          il(i+1) = il(i) + il(i+1)
          iu(i+1) = iu(i) + iu(i+1)
       END DO
       ! now do the actual copying
       DO i=1,n
          DO k=ia(i),idiag(i)-1
             j = ja(k)
             next = iu(j)
             ao(next) = a(k)
             jao(next) = i
             iu(j) = next+1
          END DO
          ao(i)=a(idiag(i))
          DO k=idiag(i)+1,ia(i+1)-1
             j = ja(k)
             next = il(j)
             ao(next) = a(k)
             jao(next) = i
             il(j) = next+1
          END DO
       END DO
       ! reshift il, iu and leave
       DO i=n,1,-1
          il(i+1) = il(i)
          iu(i+1) = iu(i)
       END DO
       il(1)=ili
       iu(1)=iui
    END IF
    RETURN
  END SUBROUTINE dagmg_csrdlu
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of CSR TO DLU FORMAT
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! MATRIX VECTOR MULTIPLICATION
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrmv(n, a, ja, ifja, ia, x, ifx, y, iad, flop, lstin)
!
!   Perform multiplication of vector x by matrix A in (n,a,ja,ia)
!   in csr format; result in y
!
!   iad==0,1 : y =  A*x
!   iad==-1  : y = -A*x
!   iad > 1  : y := y+A*x
!   iad <-1  : y := y-A*x
!
!   ja declared as ja(ifja:*) for compatibility with d+l+u format
!   (can serve for multiplication by lower or upper part setting ifja=n)
!   standard usage: ifja=1
!
!   x declared as x(ifx:*) for the need of the parallel version;
!   standard usage: ifx=1
!
!   standard usage: lstin <0
!     otherwise, lstin decalred as lstin(0:*), and
!     consider only rows listed in lstin(1:lstin(0))
!
!   flop : floating point ops counter
!
    !
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iad, ifx, i, k, j
    INTEGER, OPTIONAL :: lstin(0:*)
    REAL(kind(0.0d0)) :: a(*), x(ifx:*), y(n), t
    REAL(kind(0.0d0)) :: flop
    !
       IF (iad .LT. -1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad .GT. 1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad .EQ. -1) THEN
          DO i=1,n
             t=0.0d0
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE
          DO i=1,n
             t=0.0d0
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       END IF
       !
    flop=flop+dble(2*(ia(n+1)-ia(1)))
    RETURN
  END SUBROUTINE dagmg_csrmv
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of MATRIX VECTOR MULTIPLICATION
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! TRIANGULAR SOLVES
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrlsolve(n, a, ja, ifja, ia, p, x, y, iunit, flop)
!
!   Perform triangular solve with lower triangular matrix L;
!     srtrict lower part in (n,a,ja,ia) in csr format and
!     INVERSE of the diagonal in p(1:n)
!     result in x
!
!   x and y may point to the same memory location
!
!   iunit >0 : x = inv(L) * y
!   iunit <0 : x = inv( inv(diag(L))*L ) * y = inv(P*L) * y
!   iunit==0 : any of the above assuming unit diagonal (p not accessed)
!
!   ja declared as ja(ifja:*) for compatibility with d+l+u format
!
! flop : floating point ops counter
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    REAL(kind(0.0d0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0d0)) :: flop
!
    IF (iunit .LT. 0) THEN
       x(1) = y(1)
       DO i=2,n
          t = 0.0d0
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = y(i) + p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE IF (iunit .GT. 0) THEN
       x(1) = p(1)*y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(1) = y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE dagmg_csrlsolve
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_csrusolve(n, a, ja, ifja, ia, p, x, y, iunit,flop)
!
!   Perform triangular solve with upper triangular matrix U;
!     srtrict upper part in (n,a,ja,ia) in csr format and
!     INVERSE of the diagonal in p(1:n)
!     result in x
!
!   x and y may point to the same memory location
!
!   iunit >0 : x = inv(U) * y
!   iunit <0 : x = inv( inv(diag(U))*U ) * y = inv(P*U) * y
!   iunit==0 : any of the above assuming unit diagonal (p not accessed)
!
!   ja declared as ja(ifja:*) for compatibility with d+l+u format
!
! flop : floating point ops counter
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    REAL(kind(0.0d0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0d0)) :: flop
!
    IF (iunit .LT. 0) THEN
       x(n) = y(n)
       DO i=n-1,1,-1
          t = 0.0d0
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = y(i) + p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE IF (iunit .GT. 0) THEN
       x(n) = p(n)*y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(n) = y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE dagmg_csrusolve
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of TRIANGULAR SOLVES
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! ROUTINES FOR ITERATIVE SOLUTION
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_FlexCG(N,f,X,ITER,RESID,a,ja,ia,init,flop)
    !
    !  Implements flexible conjugate gradient iterative solver
    !   with MG preconditioner
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, init
    REAL(kind(0.0d0)) :: flop
    REAL(kind(0.0d0)) :: f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    REAL(kind(0.0d0)) :: a(*)
    INTEGER       :: MAXIT, ierr, kk
    REAL(kind(0.0d0)) :: TOL, BNORM, RESID, dum0
    REAL(kind(0.0d0)) :: ALPHA, BET0, RHO, dum(6), td
    REAL(kind(0.0d0)) :: BET1
    REAL(kind(0.0d0)), ALLOCATABLE :: SD(:)
    REAL(kind(0.0d0)), ALLOCATABLE :: S(:),fsc(:)
    REAL(kind(0.0d0)), external :: DDOT
    REAL(kind(0.0d0)), external :: DNRM2
    INTEGER , parameter :: IONE=1
    !
    mritr=nwrkcum+4*N
    ALLOCATE (S(2*N+1:nwrkcum+5*N),SD(N))
    !
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE(iout,940) IRANK
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (smoothtype == 1) THEN
          WRITE(iout,945) nsmooth
       ELSE
          WRITE(iout,946) nsmooth
       END IF
       WRITE(iout,947)
    END IF
    !
    TOL = RESID
    MAXIT = ITER
    RESID = 1.0d0
    ITER = 0
    dum(3) = DNRM2(N, f, IONE)**2
    flop=flop+dble(2*N)
    !
    !    compute initial residual
    !
    IF (init==1) THEN
       CALL dagmg_matv(N,x,SD,a,ja,ia,flop,trans )
       f(1:n)=f(1:n)-SD(1:N)
       dum(2) = DNRM2(N, f, IONE)**2
       BNORM=SQRT(dum(3))
       RESID=SQRT(dum(2))
       IF (BNORM.EQ.0.0d0) THEN
          !
          ! R.h.s. is trivial: exit with solution equal to the zero vector
          !
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=0.0d0
          RETURN
       END IF
       RESID=RESID/BNORM
       IF(wff.AND. (MAXIT <= 0 .OR. RESID <= TOL)) THEN
          WRITE(iout,900) 0, resid*bnorm, resid
       END IF
    END IF
    !
    !
    DO WHILE ( ITER < MAXIT .AND. RESID > TOL )
       ITER = ITER + 1
       !
       !        Preconditioner Solve.
       !
       CALL dagmg_prec_matv(1,f,S(1+3*N),S(1+4*N),flop,S(1+5*N),.FALSE.)
       !
       !        Compute direction vector.
       !
       IF ( ITER > 1 ) THEN
          dum(1) = - DDOT(N, SD, IONE, S(1+3*N), IONE)
          BET0=dum(1)
          BET1=BET0/RHO
          CALL DAXPY(N, BET1, S(1+2*N), IONE, S(1+3*N), IONE)
          flop=flop+dble(4*N)
       ENDIF
       CALL DCOPY(N, S(1+3*N), IONE, S(1+2*N), IONE)
       CALL dagmg_matv(N,S(1+2*N),SD,a,ja,ia,flop,trans )
       !
       dum(1) =  DDOT(N, S(1+2*N), IONE, SD, IONE)
       dum(2) =  DDOT(N,S(1+2*N),IONE,f,IONE)
       flop=flop+dble(4*N)
       !
       IF (ITER==1) THEN
          IF (init == 0) THEN
             BNORM=SQRT(dum(3))
             IF (BNORM.EQ.0.0d0) THEN
                !
                ! R.h.s. is trivial: exit with solution equal to the zero vector
                !
                IF(wff) THEN
                   WRITE(iout,998)
                   WRITE(iout,999)
                END IF
                X(1:N)=0.0d0
                RETURN
             END IF
          END IF
          IF(wff) THEN
             WRITE(iout,900) 0, resid*bnorm, resid
          END IF
       ELSE
       END IF
       !
       RHO=dum(1)
       ALPHA=dum(2)/RHO
       !
       IF (ITER == 1 .AND. init == 0) THEN
          CALL DCOPY(N,S(1+2*N),IONE,X,IONE)
          CALL DSCAL(N,ALPHA,X,IONE)
          flop=flop+dble(N)
       ELSE
          CALL DAXPY(N, ALPHA, S(1+2*N), IONE, X, IONE)
          flop=flop+dble(2*N)
       END IF
       CALL DAXPY(N, -ALPHA, SD, IONE, f, IONE)
       !
       dum0 = DNRM2(N,f,IONE)**2
       RESID=dum0
       RESID=SQRT(RESID)
       RESID=RESID/BNORM
       IF (wff) THEN
          WRITE(iout,900) iter, resid*bnorm, resid
       END IF
       flop=flop+dble(4*N)
    END DO
    !
    IF( resid > tol ) THEN
       IF (IRANK <= 0) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    !
    kstat(2,1)=ABS(iter)
    !
    DEALLOCATE(S,SD)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    RETURN
900 FORMAT('****  Iter=',i5,'        Resid=',e9.2,                 &
         '        Relat. res.=',e9.2)
940 FORMAT(i3,                                                     &
         '*SOLUTION: flexible conjugate gradient iterations (FCG(1))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
945 FORMAT(  '****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre-and post-smoothing step(s)')
946 FORMAT(  '****     AMG preconditioner with',i2,             &
         ' ILU(0) pre-and post-smoothing step(s)')
947 FORMAT(  '****         at each level')
950 FORMAT('**** !!!   WARNING!!!',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
    !
  END SUBROUTINE dagmg_FlexCG
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_GCR(N,f,X,ITER,RESID,a,ja,ia,init,flop)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, init
    REAL(kind(0.0d0)) ::  flop
    REAL(kind(0.0d0)) ::  f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    REAL(kind(0.0d0)) :: a(*)
    INTEGER   :: MAXIT,i,itm,irst,ierr,j,k
    REAL(kind(0.0d0)) ::  ALPHA, BET0
    REAL(kind(0.0d0)) ::  RESID, RESID2, BNORM2,  RHO, TRS
    REAL(kind(0.0d0)) ::  TOL, TOL2BNORM2, TOLT, dum0
    REAL(kind(0.0d0)) :: dum(6),xd,t
    REAL(kind(0.0d0)), ALLOCATABLE :: Sc(:),Siv(:),SiR(:)
    REAL(kind(0.0d0)), ALLOCATABLE :: Su(:),R(:),fsc(:),Sw(:)
    REAL(kind(0.0d0)), external :: DDOT
    REAL(kind(0.0d0)), external :: DNRM2
    INTEGER , parameter :: IONE=1
    INTEGER  :: itm1, m, info
    !
    mritr=nwrkcum+N*2*nrst+((nrst+1)*nrst)/2+nrst
    ALLOCATE (Su(N*nrst),Sc(N*nrst)                  &
             ,SiR(((nrst+1)*nrst)/2),Siv(nrst)       &
             ,R(nwrkcum))
    !
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE(iout,938) IRANK,nrst
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (smoothtype == 1) THEN
          WRITE(iout,945) nsmooth
       ELSE
          WRITE(iout,946) nsmooth
       END IF
       WRITE(iout,947)
    END IF
    !
    m=MIN(nrst,ITER)
    TRS=EPSILON(1.0d0)
    TRS=SQRT(SQRT(TRS))
    TOL = RESID
    TOL2BNORM2 = TOL
    MAXIT = ITER
    RESID2= 1.0d0
    ITER = 0
    dum(3) = DNRM2(N, f, IONE)**2
    flop=flop+dble(2*N)
    !
    !    compute initial residual
    !
    IF (init==1) THEN
       CALL dagmg_matv(N,x,Sc,a,ja,ia,flop,trans )
       f(1:n)=f(1:n)-Sc(1:N)
       dum(2) = DNRM2(N, f, IONE)**2
       BNORM2=dum(3)
       RESID2=dum(2)
       IF (BNORM2.EQ.0.0d0) THEN
          !
          ! R.h.s. is trivial: exit with solution equal to the zero vector
          !
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=0.0d0
          RETURN
       END IF
       TOL2BNORM2 = TOL*TOL*BNORM2
       IF (wff.AND. (MAXIT <= 0 .OR. RESID2 <= TOL2BNORM2)) THEN
          WRITE(iout,900) 0, 0,SQRT(resid2),SQRT(resid2/bnorm2)
       END IF
    END IF
    !
    itm  = -1
    irst = 0
    DO WHILE ( ITER < MAXIT .AND. RESID2 > TOL2BNORM2 )
       itm  = itm  + 1
       ITER = ITER + 1
       !
       !      Restarting
       IF (itm == m) THEN
          CALL DTPTRS('U','N','U',m,IONE,SiR,Siv,m,info)
          IF (irst == 0 .AND. init == 0) THEN
             CALL DGEMV('N',N,m,1.0d0,Su,       &
                  N,Siv,IONE,0.0d0,X,IONE)
             flop=flop+dble(2*m*N+m*(m+1))
          ELSE
             CALL DGEMV('N',N,m,1.0d0,Su,        &
                  N,Siv,IONE,1.0d0,X,IONE)
             flop=flop+dble((2*m+1)*N+m*(m+1))
          END IF
          itm=0
          irst=irst+1
       END IF
       !
       !        Preconditioner Solve & MATVEC
       !
       CALL dagmg_prec_matv(1,f,Su(1+itm*N),Sc(1+itm*N),flop,R,.FALSE.)
       CALL dagmg_matv(N, Su(1+itm*N), Sc(1+itm*N)             &
            , a, ja, ia, flop,trans )
       !
       !        Gram-Schmidt
       !
       IF (itm > 0) THEN
          DO i=0,itm-1
             dum(1)=DDOT(N,Sc(1+i*N),IONE,Sc(1+itm*N),IONE)
             bet0=dum(1)
             bet0=bet0/SiR(1+i+(i*(i+1))/2)
             SiR(1+i+(itm*(itm+1))/2)=bet0
             CALL DAXPY(N,-bet0,Sc(1+i*N),IONE,Sc(1+itm*N),IONE)
             flop=flop+dble(4*N)
          END DO
       END IF
       !
       !           no normalisation: record norm instead
       !
       dum(1)=DNRM2(N,Sc(1+itm*N),IONE)**2
       dum(2)=DDOT(N, Sc(1+itm*N), IONE, f, IONE)
       IF (ITER == 1) THEN
          IF (init == 0) THEN
             BNORM2=dum(3)
             RESID2=BNORM2
             IF (BNORM2.EQ.0.0d0) THEN
                !
                ! R.h.s. is trivial: exit with solution equal to the zero vector
                !
                IF(wff) THEN
                   WRITE(iout,998)
                   WRITE(iout,999)
                END IF
                X(1:N)=0.0d0
                RETURN
             END IF
             TOL2BNORM2=TOL*TOL*BNORM2
          END IF
          IF (wff) THEN
             WRITE(iout,900) 0, 0,SQRT(resid2),SQRT(resid2/bnorm2)
          END IF
          TOLT = MAX(TOL2BNORM2,TRS*RESID2)
       ELSE
       END IF
       rho=dum(1)
       alpha=dum(2)
       SiR(1+itm+(itm*(itm+1))/2)=rho
       bet0=alpha/rho
       Siv(1+itm)=bet0
       !
       CALL DAXPY(N, -bet0, Sc(1+itm*N), IONE, f, IONE)
       flop=flop+dble(6*N)
       !
       RESID2 = RESID2 - alpha*alpha/rho
       IF (RESID2 <= TOLT) THEN
          RESID2 = DNRM2(N,f,IONE)**2
          flop=flop+dble(2*N)
          TOLT = MAX(TOL2BNORM2,TRS*RESID2)
       END IF
       IF (wff)THEN
          WRITE(iout,900) iter,irst,SQRT(ABS(resid2)),SQRT(ABS(resid2/bnorm2))
       END IF
       !
    END DO
    !
    IF (itm >= 0) THEN
       itm1=itm+1
       CALL DTPTRS('U','N','U',itm1, IONE,SiR,Siv,m,info)
       IF (irst == 0 .AND. init == 0) THEN
          CALL DGEMV('N',N,itm1,1.0d0,Su,        &
               N,Siv,IONE,0.0d0,X,IONE)
          flop=flop+dble(2*(itm+1)*N+(itm+1)*(itm+2))
       ELSE
          CALL DGEMV('N',N,itm1,1.0d0,Su,        &
               N,Siv,IONE,1.0d0,X,IONE)
          flop=flop+dble((2*(itm+1)+1)*N+(itm+1)*(itm+2))
       END IF
    END IF
    !
    RESID=SQRT(RESID2/BNORM2)
    IF( resid > tol ) THEN
       IF (IRANK <= 0) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    !
    DEALLOCATE (Su,Sc,R,Siv,SiR)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    IF(ALLOCATED(Sw)) DEALLOCATE(Sw)
    !
    kstat(2,1)=ABS(iter)
    !
    RETURN
900 FORMAT('****  Iter=',i5,' (',i2,' rest.)        Resid=',e9.2,    &
         '        Relat. res.=',e9.2)
938 FORMAT(i3,'*SOLUTION: GCR iterations (GCR(',i2,'))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
945 FORMAT(  '****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre-and post-smoothing step(s)')
946 FORMAT(  '****     AMG preconditioner with',i2,             &
         ' ILU(0) pre-and post-smoothing step(s)')
947 FORMAT(  '****         at each level')
950 FORMAT('**** !!!   WARNING!!!',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
    !
  END SUBROUTINE dagmg_GCR
!--------------------------------------------------------------------
  SUBROUTINE dagmg_matv(n, x, y, a, ja, ia, flop, transpose,         &
       iext, lstout, ilstout, lstin, ilstin)
!
!  Performs matrix vector product by x, result in y
!  If trans==.TRUE., multiply by the transpose
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), i, kk, k1, k2, ier, idum(1)
    INTEGER, OPTIONAL :: iext(n), lstout(*), ilstout(*), lstin(0:*), ilstin(*)
    REAL(kind(0.0d0)) :: y(n), a(*), t, xx
    REAL(kind(0.0d0)) :: x(n)
    REAL(kind(0.0d0)) :: flop
    LOGICAL :: transpose
    !
    !
    IF (.NOT.transpose) THEN
       DO i = 1, n
          k1 = ia(i)
          xx = x(ja(k1))
          t = a(k1) * xx
          k2 = ia(i+1)-1
          DO kk = k1+1, k2
             xx = x(ja(kk))
             t = t + a(kk)*xx
          ENDDO
          y(i) = t
       ENDDO
    ELSE
       y(1:n)=0.0d0
       DO i = 1, n
          xx = x(i)
          DO kk = ia(i), ia(i+1)-1
             y(ja(kk)) = y(ja(kk)) + a(kk)*xx
          ENDDO
       ENDDO
    END IF
    !
    flop = flop + dble(2 *(ia(n+1)-1)-n)
    !
    RETURN
  END SUBROUTINE dagmg_matv
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_prec_matv(l,B,X,AX,flop,R,matv)
!
! Applies AMG preconditioner at level l to vector B ans store the result in X
! If matv==.TRUE. : computes in addition AX = A*X
! R : workspace
!     length needed at level l: lr(l)
!      nsmooth*smoothtype == 1:
!                | nn(l+1)              if l+1==nlev
!        lr(l) = | 3*nn(l+1)+lr(l+1)    otherwise, if innermax(l+1)<=1
!                | 5*nn(l+1)+lr(l+1)    otherwise
!
!      otherwise:
!                | max( 2*n , nn(l+1) )              if l+1==nlev
!        lr(l) = | max( 2*n , 3*nn(l+1)+lr(l+1) )    otherwise,
!                |                                        if innermax(l+1)<=1
!                | max( 2*n , 5*nn(l+1)+lr(l+1) )    otherwise
!
! flop : floating point ops counter
!
! AX(1:nn(l)) needed as workspace
!       (i.e., to be allocated by calling program even if matv==.FALSE.)
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l
    REAL(kind(0.0d0)), OPTIONAL :: flop
    REAL(kind(0.0d0)), OPTIONAL ::  B(*), X(*), AX(*), R(*)
    REAL(kind(0.0d0)) ::  dum(1)
    LOGICAL, OPTIONAL :: matv
    LOGICAL :: update
    INTEGER :: is,n,nnext,XN,BN,AXN,RN,idum(1)
!
!  With ILU(0) (smoothtype==2), i.e. smoothing steps implies
!  1 forward + 1 backward solve; i.e., two times more than with
!  Gauss-Seidel (smoothtype==1)
!
    INTEGER, PARAMETER :: nsmotot=nsmooth*smoothtype
    n=nn(l)
    nnext=nn(l+1)
    !
    IF (l == nlev) THEN
       !
       ! Special case: this level is the last one
       ! (Normally, one enters this section iff nlev==1)
       !
       X(1:N)=B(1:N)
       CALL dagmg_MUMPSseq(n,dt(l)%a,dt(l)%ja,dt(l)%ia,X,2,flop)
       !
       RETURN
       !
    END IF
    !
    !
    ! Presmoothing
    !
    IF (nsmotot == 1) THEN
       !
       ! Perform forward GS sweep to solve A*X= B; new residual in AX
       !  (X is new: update=.FALSE. ; R will not be used)
       CALL dagmg_fwGS(l,B,X,AX,flop,.FALSE.,R)
       !
    ELSE
       !
       ! Perform nsmotot/2 forward+backward sweeps to solve A*X=B
       !  with initial approx. in X and residual in AX when update=.TRUE.
       !  (i.e., symmetric GS or ILU(0), depending on the pivots);
       !  result in X, new residual in AX
       !  (update=.FALSE. at first call, then .TRUE.)
       ! R(1:2*n) as workspace
       update=.FALSE.
       DO is=2,nsmotot,2
          CALL dagmg_fwbwsolve1(l,dt(l)%a,B,X,AX,flop,update,R)
          update=.TRUE.
       END DO
       IF (mod(nsmotot,2) .EQ. 1) THEN
          !
          ! Perform final forward sweep to solve A*X=B
          !  with initial approx. in X and residual in AX
          !  (X has to be updated: update=.TRUE.)
          ! R(1:n) as workspace
          CALL dagmg_fwGS(l,B,X,AX,flop,.TRUE.,R)
          !
       END IF
    END IF
    !
    ! Coarse grid correction
    !
    IF (nnext > 0) THEN
       !
       ! Standard case: coarse grid is nonmpty
       !  (otherwise: nothing to be done for coarse grid correction)
       ! Restricted residual will be in R(nnext+1:2*nnext) except at last level
       ! and solution of coarse system in R(1:nnext)
       !
       XN=1
       BN=XN+nnext
       IF (l+1 == nlev) BN=XN
       ! Restric the residual in AX
       CALL dagmg_restag(N,nnext,AX,R(BN),dt(l)%ind,flop)
       !
       IF (l+1 == nlev) THEN
          !
          ! Direct solution at last level
          CALL dagmg_MUMPSseq(nnext,dt(l+1)%a,dt(l+1)%ja,dt(l+1)%ia,R,2,flop)
          !
       ELSE IF ( innermax(l+1) <= 1 ) THEN
          !
          ! V cycle enforced: apply MG preconditiner
          !   size of workspace in R at level l+1: length(R)-3*nnext
          AXN=BN+nnext
          RN=AXN+nnext
          CALL dagmg_prec_matv(l+1,R(BN),R(XN),R(AXN),flop,R(RN),.FALSE.)
          !
          kstat(1,l+1)=MAX(kstat(1,l+1),1)
          kstat(2,l+1)=kstat(2,l+1)+1
          !
       ELSE
          !
          ! K-cycle with Galerkin projection
          !  will call prec_matv at level l+1 with workspace in R of size
          !  at least length(R(BN:*))-4*nnext, i.e. length(R)-5*nnext
          CALL dagmg_inner_iter(nnext,R(XN),R(BN),l+1,flop)
          !
       END IF
       !
       ! Prolongate just computed coarse solution, and add it to X
       CALL dagmg_prolag(N,nnext,X,R(XN),dt(l)%ind,flop)
       !
    END IF
    !
    ! Post smoothing
    !
    IF (nsmotot == 1) THEN
       !
       ! Perform backward GS sweep to solve B-A*X  with initial approx. in X
       !  AX to be computed if matv==.TRUE.
       CALL dagmg_bwGS(l,B,X,AX,flop,matv)
       !
    ELSE
       IF (mod(nsmotot,2) .EQ. 1) THEN
          !
          !  Perform initial backward sweep to solve B-A*X
          !   with initial approx. in X; do not compute AX
          !   (but use it as worskspace)
          CALL dagmg_bwGS(l,B,X,AX,flop,.FALSE.)
          !
       END IF
       !
       ! Apply nsmotot/2 forward+backward sweeps to solve A*X=B
       !  with initial approx. in X
       !  (i.e., symmetric GS or ILU(0), depending on the pivots)
       !  compute AX at last call if required
       ! R(1:2*n) as workspace
       update=.FALSE.
       DO is=2,nsmotot,2
          IF (is .GE. nsmotot-1) update=matv
          CALL dagmg_fwbwsolve2(l,dt(l)%a,B,X,AX,flop,update,R)
       END DO
    END IF
    RETURN
  END SUBROUTINE dagmg_prec_matv
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_fwGS(l,B,X,AX,flop,update,R)
!
! Performs pre-smoothing step at level l with smoother:
!
!          M = (D+E)
!
!          E strictly lower triangular,
!          D diagonal
!
!   Elements of E stored in (dt(l)%a,dt(l)%ja,dt(l)%il)
!   Inverse of the elements of D strored in dt(l)%a(1:nn(l))
!
!   If update==.FALSE., computes X = inv(M)*B and AX = B - A*X,
!   based on the assumption that either A=D+E+F with elements of F
!   stored in (dt(l)%a,dt(l)%ja,dt(l)%iu) (strict upper part)
!
!   If update==.TRUE., computes X := X + inv(M)*AX and AX := AX - A-inv(M)*AX,
!   based on the assumption that either A=D+E+F with elements of F
!   stored in (dt(l)%a,dt(l)%ja,dt(l)%iu) (strict upper part)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) :: flop
    REAL(kind(0.0d0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: update
    n=nn(l)
    !
    IF (update) THEN
       !
       ! Solve (D+E) R = AX
       !
       CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a   &
                            ,R,AX,1,flop)
       X(1:n)=X(1:n)+R(1:n)
       flop=flop+dble(n)
       !
       ! AX_old - A*R = AX_old - F * R - (D+E+) * inv(D+E) * AX_old = - F * R
       !
       ! Hence, compute AX = - F * R
       !
       CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,R,1,AX,-1,flop)
       !
    ELSE
       !
       ! Solve (D+E) X = B
       !
       CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a,X,B,1,flop)
       !
       ! B - A*X = B - F * X - (D+E) * inv(D+E) * B = - F * X
       ! Hence, compute AX = - F * R
       !
       CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,-1,flop)
       !
    END IF
    RETURN
  END SUBROUTINE dagmg_fwGS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_bwGS(l,B,X,AX,flop,matv)
!
! Performs post-smoothing step at level l with smoother:
!
!          M = (D+F)
!
!          F strictly upper triangular,
!          D diagonal
!
!   Elements of E stored in (dt(l)%a,dt(l)%ja,dt(l)%iu)
!   Inverse of the elements of D strored in dt(l)%a(1:nn(l))
!
!   Computes X := X + inv(M)*(B-A*X),
!   based on the assumption that either A=D+E+F with elements of E
!   stored in (dt(l)%a,dt(l)%ja,dt(l)%il) (strict lower part)
!
!   If matv==.TRUE., computes in addition AX=A*X,
!   based on the same assumptions
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) :: flop
    REAL(kind(0.0d0)) ::  B(*), X(*), AX(*)
    LOGICAL :: matv
    n=nn(l)
    !
    ! Compute AX = B - E * X
    !
    AX(1:n)=B(1:n)
    CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
    !
    ! Solve (D+F) X = AX
    ! Then:  X := inv((D+F)) * ( B - E * X ) = X + inv(D+F) * ( B - A * X )
    !
    CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a,X,AX,1,flop)
    !
    IF (.NOT.matv) RETURN
    !
    ! Computation of AX := A*X is needed
    ! Use A*X = (D+F) * X + E * X = AX + E * X
    !
    CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    RETURN
  END SUBROUTINE dagmg_bwGS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_fwbwsolve1(l,p,B,X,AX,flop,update,R)
!
! Performs pre-smoothing step at level l with smoother:
!
!          M = (P+E) inv(P) (P+F)
!
!          E strictly lower triangular,
!          F strictly upper triangular,
!          P diagonal
!
!   Elements of E stored in (dt(l)%a,dt(l)%ja,dt(l)%il)
!   Elements of F stored in (dt(l)%a,dt(l)%ja,dt(l)%iu)
!   Inverse of the elements of P strored in p(1:nn(l))
!
!   If update==.FALSE., computes X = inv(M)*B and AX = B - A*X,
!   based on the assumption that either (smoothtype==1) A=P+E+F
!   or (smoothtype==2) A=P+E+F+Q with Q diagonal stored in
!   dt(l)%a(1:n)
!
!   If update==.TRUE., computes X := X + inv(M)*AX and AX := AX - A-inv(M)*AX,
!   based on the assumption that either (smoothtype==1) A=P+E+F
!   or (smoothtype==2) A=P+E+F+Q with Q diagonal stored in
!   dt(l)%a(1:n)
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) ::  p(*), B(*), X(*), AX(*), R(*)
    REAL(kind(0.0d0)) :: flop
    LOGICAL :: update
    !
    n=nn(l)
    IF (.NOT.update) THEN
       !
       ! Solve (I+E*inv(P)) R(1:n) = B
       !
       CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,p    &
                            ,R,B,-1,flop)
       !
       ! Solve (P+F) X = R(1:n)
       !
       CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,p    &
                             ,X,R,1,flop)
       !
       ! B - A*X = B - (P+E+F+Q) * inv(P+F) * R(1:n)
       !         = B - R(1:n) - (E+Q)*X
       !
       !
       ! First compute AX = B - R(1:n)
       !
       AX(1:n)=B(1:n)-R(1:n)
       flop=flop+dble(n)
       !
       ! Then AX := AX - E*X
       !
       CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
       !
       ! Eventually, AX := AX - Q*X if needed
       IF (smoothtype == 2) AX(1:n)=AX(1:n)-dt(l)%a(1:n)*X(1:n)
       !
    ELSE
       !
       ! Solve (I+E*inv(P)) R = AX
       !
       CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,p    &
                            ,R,AX,-1,flop)
       !
       ! Solve (P+F) R(n+1:2*n) = R(1:n)
       !
       CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,p    &
                            ,R(n+1),R,1,flop)
       !
       ! X = X + inv(M) AX = X + R(n+1:2*n)
       !
       X(1:n)=X(1:n)+R(n+1:2*n)
       !
       ! AX - A*inv(M)*AX = AX - (P+E+F+Q) * inv(P+F) * R(1:n)
       !                  = AX - R(1:n) - (E+Q)*R(n+1:2*n)
       !
       !
       ! First compute AX = AX - R(1:n)
       !
       AX(1:n)=AX(1:n)-R(1:n)
       flop=flop+dble(2*n)
       !
       ! Then AX := AX - E*R(n+1:2*n)
       !
       CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,R(n+1),1,AX  &
                         ,-2,flop)
       !
       ! Eventually, AX := AX - Q*R(n+1:2*n) if needed
       IF (smoothtype == 2) AX(1:n)=AX(1:n)-dt(l)%a(1:n)*R(n+1:2*n)
       IF (smoothtype == 2) flop=flop+dble(2*n)
       !
    END IF
    !
    RETURN
  END SUBROUTINE dagmg_fwbwsolve1
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_fwbwsolve2(l,p,B,X,AX,flop,matv,R)
!
! Performs post-smoothing step at level l with smoother:
!
!          M = (P+E) inv(P) (P+F)
!
!          E strictly lower triangular,
!          F strictly upper triangular,
!          P diagonal
!
!   Elements of E stored in (dt(l)%a,dt(l)%ja,dt(l)%il)
!   Elements of F stored in (dt(l)%a,dt(l)%ja,dt(l)%iu)
!   Inverse of the elements of P strored in p(1:nn(l))
!
!   Computes X := X + inv(M)*(B-A*X),
!   based on the assumption that either (smoothtype==1) A=P+E+F
!   or (smoothtype==2) A=P+E+F+Q with Q diagonal stored in
!   dt(l)%a(1:n)
!
!   If matv==.TRUE., computes in addition AX=A*X,
!   based on the same assumptions
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) ::  p(*), B(*), X(*), AX(*), R(*)
    REAL(kind(0.0d0)) :: flop
    LOGICAL :: matv
    !
    n=nn(l)
    !
    ! inv(M) * (B-A*X) = inv(I+inv(P)*F) * inv(P+E) *( B- (P+E+F+Q)*X )
    !                  = inv(I+inv(P)*F) * ( inv(P+E) * (B-(F+Q)*X) - X )
    !
    ! First compute AX=F*X
    !
    CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,0,flop)
    !
    ! Then, R(1:n)= B - AX
    R(1:n)=B(1:n)-AX(1:n)
    !    and R(1:n) := R(1:n) - Q*X if needed
    IF (smoothtype == 2) R(1:n)=R(1:n)-dt(l)%a(1:n)*X(1:n)
    IF (smoothtype == 2) flop=flop+dble(2*n)
    !
    !
    ! Now solve (P+E) R(n+1:2*n)=R(1:n)
    !
    CALL dagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,p    &
                         ,R(n+1),R,1,flop)
    !
    ! and solve (I+inv(P)*F) R(1:n) = R(n+1:2*n) - X
    !
    R(1:n) = R(n+1:2*n) - X(1:n)
    CALL dagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,p    &
                         ,R,R,-1,flop)
    !
    ! Eventually, X := X + inv(M) * (B-A*X) = X + R(1:n)
    !
    X(1:n)=X(1:n) + R(1:n)
    flop=flop+dble(3*n)
    !
    IF (.NOT.matv) RETURN
    ! If needed, compute
    !      AX := A * X
    !          = (P+E+F+Q) * ( X_old + R(1:n) )
    !          = AX + (P+E+Q)* X_old + (E+Q) * R(1:n)
    !                 (P+F) * inv((I+inv(P)*F)) * (R(n+1:2*n)-X_old)
    !          = AX + (P+E+Q)*X_old + (E+Q) * R(1:n) + P * (R(n+1:2*n)-X_old)
    !          = AX + P * R(n+1:2*n) + (E+Q) * X
    !
    AX(1:n)=AX(1:n)+R(n+1:2*n)/P(1:n)
    flop=flop+dble(2*n)
    !
    ! Now, AX: = AX + E*X
    !
    CALL  dagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    !     and AX := AX + Q*X if needed
    IF (smoothtype == 2) AX(1:n)=AX(1:n)+dt(l)%a(1:n)*R(n+1:2*n)
    IF (smoothtype == 2) flop=flop+dble(2*n)
    RETURN
  END SUBROUTINE dagmg_fwbwsolve2
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_inner_iter(n,X,R,l,flop)
    !
    ! Perform inner iterations (K-cycle) for AGMG
    ! Based on projection:
    !   Gelerkin projection if gcrin==.FALSE.; in the SPD case, this is
    !                  equivalent to Flexible CG, implemented as in [1];
    !   Minimal residual projection if gcrin==.TRUE.; i.e.,
    !                  GCR implemented as in [1].
    !   In addition, peforms automatic switch to Minimal residual projection
    !           if the projection of symmetric part of the matrix is
    !           checked to be not postive definite.
    !
    ! Input: N       : number of unknowns at that level
    !        R(1:N)  : r.h.s. of the system  (modified on output)
    !        R(N+1,*): work space
    !              needed length for R: 4*n + workspace for prec_matv
    !
    !              [ comments below: Ri means R(1:n,i) ]
    !
    !        l       : level index
    !
    ! Output: X(1:n) : computed solution
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER   :: N, ITER, l, ierr, k
    REAL(kind(0.0d0)) ::  RESID, flop, BNORM, det
    REAL(kind(0.0d0)) :: X(N), R(N,*)
    REAL(kind(0.0d0)) :: alpha1,alpha2,bet0,bet1,bet2,rho1,rho2,gamm0,gamm1,zeta
    REAL(kind(0.0d0)) :: dum(10),y1,y2
    REAL(kind(0.0d0)), external :: DDOT
    REAL(kind(0.0d0)), external :: DNRM2
    INTEGER , parameter :: IONE=1
    !
    ! AT MOST 2 ITERATIONS
    !
    dum(3)=DNRM2(N, R, IONE)**2
    !
    ! Early exit if r.h.s. is trivial
    !
    IF (dum(3) .EQ. 0.0d0) THEN
       X(1:N)=0.0d0
       flop=flop+dble(2*N)
       RETURN
    END IF
    ITER = 1
    !
    !  Apply preconditioner to r.h.s. (result in X) and compute
    !    matrix vector pruduct (result in R2: R2=A*X)
    !    workspace tranferred to prec_matv: length(R)-2*n
    !
    CALL dagmg_prec_matv(l,R,X,R(1,2),flop,R(1,3),.TRUE.)
    !
!!$    IF (gcrin) GOTO 100
    !
    dum(1) = DDOT(N,X,IONE,R(1,2),IONE)
    dum(2) = DDOT(N,X,IONE,R,IONE)
    IF (ABS(dum(1)) .LE. repsmach*ABS(dum(2))) THEN
       ! Projected matrix near singular: switch to GCR projection
!!$       gcrin=.TRUE.
       flop=flop+dble(2*N)
       GOTO 100
    END IF
    BNORM = dum(3)
    rho1=dum(1)
    alpha1=dum(2)
    !
    bet0=alpha1/rho1
    !
    !
    !  Compute residual corresponding to solution bet0 * X
    CALL DAXPY(N, -bet0, R(1,2), IONE, R, IONE)
    !  ... and the corresponding square norm
    RESID = DNRM2(N,R,IONE)**2
    IF (RESID <= resi*resi*BNORM) THEN
       ! criterion satified, exit with solution bet0 * X
       ! Galerkin projection: the solution computed in this way
       !   is such that the corresponding residual
       !      R - bet0*A*X  = R - bet0*R2
       !   is orthogonal to X
       !   Indeed: bet0= ( X , R ) / ( X , R2 )
       CALL DSCAL( N, bet0, X, IONE )
       flop=flop+dble(11*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    !
    ITER = 2
    !
    !  Apply preconditioner to this new residual (result in R3) and compute
    !    matrix vector pruduct (result in R4: R4=A*R3)
    !  R3 is thus the second direction vector
    !    workspace tranferred to prec_matv: length(R)-4*n
    !
    CALL dagmg_prec_matv(l,R,R(1,3),R(1,4),flop,R(1,5),.TRUE.)
    !
    !
!!$    ! Check if gcrin has been switched while performing iter at lower levels
!!$    IF (gcrin) GOTO 200
    !
    dum(1) = DDOT(N,R(1,3),IONE,R(1,2),IONE)
    dum(2) = DDOT(N,R(1,3),IONE,R,IONE)
    dum(3) = DDOT(N,R(1,3),IONE,R(1,4),IONE)
    IF (.NOT.spd) dum(4) = DDOT(N,X,IONE,R(1,4),IONE)
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    IF (spd) THEN
       ! in the SPD case:
       ! gamm1 := ( X , R4 ) = ( X, A*R3 ) = ( A*X , R3 ) = ( R2 , R3)
       gamm1 = gamm0
    ELSE
       gamm1 = dum(4)
    END IF
    !
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    !
    IF (ABS(bet1).LE.repsmach*ABS(alpha2) .OR.   &
        ABS(bet2)*repsmach.GE.1.0d0)            THEN
       ! Projected matrix near singular: switch to GCR projection
!!$       gcrin=.TRUE.
       flop=flop+dble(6*N)
       IF (.NOT.spd) flop=flop+dble(2*N)
      GOTO 200
    END IF
    CALL DSCAL(N, bet2, X, IONE)
    CALL DAXPY(N, zeta, R(1,3), IONE, X, IONE)
    !
    ! Galerkin projection: the solution computed in this way
    !   is such that the corresponding residual
    !      R - bet2*A*X - zeta*A*R3 = R - bet2*R2 - zeta*R4
    !   is orthogonal to X and R3
    !   Indeed:
    !       ( X , R ) - bet2* ( X , R2 ) - zeta* ( X , R4 )
    !          =  alpha1 - bet2*rho1 - zeta*gamm1
    !          =  alpha1 - (alpha1-alpha2*gamm1/bet1) - (alpha2/bet1)*gamm1
    !          =  0
    !       ( R3 , R ) - bet2* ( R3 , R2 ) - zeta* ( R3 , R4 )
    !           [!: alpha2 is the dot product with update residual ]
    !          = ( R3 , R-bet0*R2) + bet0* ( R3 , R2)
    !                              - bet2* ( R3 , R2 ) - zeta* ( R3 , R4 )
    !          =  alpha2 + bet0*gamm0 - bet2*gamm0 - zeta*rho2
    !          =  alpha2 + alpha1*gamm0/rho1
    !                    - (alpha1-alpha2*gamm1/bet1)*(gamm0/rho1)
    !                    - (alpha2/bet1)*(bet1+gamm0*gamm1/rho1)
    !          =  0
    !
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(19*N)
    IF (.NOT.spd) flop=flop+dble(2*N)
    !
    RETURN
    !
100 CONTINUE
    dum(1)=DNRM2(N,R(1,2),IONE)**2
    dum(2)=DDOT(N, R(1,2), IONE, R, IONE )
    BNORM = dum(3)
110 CONTINUE
    rho1=dum(1)
    alpha1=dum(2)
    !
    bet0=alpha1/rho1
    !
    ! RESID below is the square of the residual norm for solution bet0 * X
    RESID=BNORM-alpha1*bet0
    IF (RESID <= resi*resi*BNORM) THEN
       ! criterion satified, exit with solution bet0 * X
       ! GCR projection: the solution computed in this way
       !   is such that the corresponding residual
       !      R - bet0*A*X  = R - bet0*R2
       !   is orthogonal to R2
       !   Indeed: bet0= ( R2 , R ) / ( R2 , R2 )
       CALL DSCAL( N, bet0, X, IONE )
       flop=flop+dble(7*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    !
    !  Compute residual corresponding to solution bet0 * X
    CALL DAXPY(N, -bet0, R(1,2), IONE, R, IONE)
    !
    ITER = 2
    !
    !  Apply preconditioner to this new residual (result in R3) and compute
    !    matrix vector pruduct (result in R4: R4=A*R3)
    !  R3 is thus the second direction vector
    !    workspace tranferred to prec_matv: length(R)-4*n
    !
    CALL dagmg_prec_matv(l,R,R(1,3),R(1,4),flop,R(1,5),.TRUE.)
    !
    !
    dum(1) = DDOT(N,R(1,4),IONE,R(1,2),IONE)
    dum(2) = DDOT(N,R(1,4),IONE,R,IONE)
    dum(3) = DNRM2(N,R(1,4),IONE)**2
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    ! GCR case:
    ! gamm1 := ( R2 , R4 ) = ( A*X, A*R3 ) = ( A*R3 , A*X ) = ( R4 , R2)
    gamm1 = gamm0
    !
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    !
    CALL DSCAL(N, bet2, X, IONE)
    CALL DAXPY(N, zeta, R(1,3), IONE, X, IONE)
    !
    ! GCR projection: the solution computed in this way
    !   is such that the corresponding residual
    !      R - bet2*A*X - zeta*A*R3 = R - bet2*R2 - zeta*R4
    !   is orthogonal to R2 and R4
    !   Indeed:
    !       ( R2 , R ) - bet2* ( R2 , R2 ) - zeta* ( R2 , R4 )
    !          =  alpha1 - bet2*rho1 - zeta*gamm1
    !          =  alpha1 - (alpha1-alpha2*gamm1/bet1) - (alpha2/bet1)*gamm1
    !          =  0
    !       ( R4 , R ) - bet2* ( R4 , R2 ) - zeta* ( R4 , R4 )
    !           [!: alpha2 is the dot product with update residual ]
    !          = ( R4 , R-bet0*R2) + bet0* ( R4 , R2)
    !                              - bet2* ( R4 , R2 ) - zeta* ( R4 , R4 )
    !          =  alpha2 + bet0*gamm0 - bet2*gamm0 - zeta*rho2
    !          =  alpha2 + alpha1*gamm0/rho1
    !                    - (alpha1-alpha2*gamm1/bet1)*(gamm0/rho1)
    !                    - (alpha2/bet1)*(bet1+gamm0*gamm1/rho1)
    !          =  0
    !
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(19*N)
    !
    RETURN
    !
200 CONTINUE
    !
    ! Switching after 1 iter: compute explicitly projected system
    !
    !   (A*X A*R2)^T * A * (X R2) * y = (A*X A*R2)^T * R
    !
    ! i.e.:
    !
    !  (R2 R4)^T * (R2 R4) * y = (R2 R4)^T * [ (R-bet0*R2) + bet0*R2 ]
    !
    ! where R-bet0*R2 is the update residual.
    !
    !  Hence
    !
    !  | dum1  dum2 |  | y1 |   | dum4 + bet0*dum1 |
    !  |            |  |    | = |                  |
    !  | dum2  dum3 |  | y2 |   | dum5 + bet0*dum2 |
    !
    !  where:
    !
    dum(1) = DNRM2(N,R(1,2),IONE)**2
    dum(2) = DDOT(N,R(1,4),IONE,R(1,2),IONE)
    dum(3) = DNRM2(N,R(1,4),IONE)**2
    dum(4) = DDOT(N,R(1,2),IONE,R,IONE)
    dum(5) = DDOT(N,R(1,4),IONE,R,IONE)
     !
     dum(4) = dum(4)+bet0*dum(1)
     dum(5) = dum(5)+bet0*dum(2)
     det = dum(1)*dum(3)-dum(2)*dum(2)
     y1  = (dum(3)*dum(4)-dum(2)*dum(5))/det
     y2  = (-dum(2)*dum(4)+dum(1)*dum(5))/det
     CALL DSCAL( N, y1, X, IONE )
     CALL DAXPY( N, y2, R(1,3), IONE, X, IONE )
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(23*N)
    !
    RETURN
    !
  END SUBROUTINE dagmg_inner_iter
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_Galerkin_inner(n,X,R,l,flop)
    !
    ! Perform inner iterations (K-cycle) for AGMG
    ! Based on Galerkin projection; in the SPD case, this is
    !   equivalent to Flexible CG, implemented as in [1]
    !
    ! Input: N       : number of unknowns at that level
    !        R(1:N)  : r.h.s. of the system  (modified on output)
    !        R(N+1,*): work space
    !              needed length for R: 4*n + workspace for prec_matv
    !
    !              [ comments below: Ri means R(1:n,i) ]
    !
    !        l       : level index
    !
    ! Output: X(1:n) : computed solution
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER   :: N, ITER, l, ierr, k
    REAL(kind(0.0d0)) ::  RESID, flop, BNORM
    REAL(kind(0.0d0)) :: X(N), R(N,*)
    REAL(kind(0.0d0)) :: alpha1,alpha2,bet0,bet1,bet2,rho1,rho2,gamm0,gamm1,zeta
    REAL(kind(0.0d0)) :: dum(10)
    REAL(kind(0.0d0)), external :: DDOT
    REAL(kind(0.0d0)), external :: DNRM2
    INTEGER , parameter :: IONE=1
    !
    ! AT MOST 2 ITERATIONS
    !
    dum(3)=DNRM2(N, R, IONE)**2
    !
    ! Early exit if r.h.s. is trivial
    !
    IF (dum(3) .EQ. 0.0d0) THEN
       X(1:N)=0.0d0
       flop=flop+dble(2*N)
       RETURN
    END IF
    ITER = 1
    !
    !  Apply preconditioner to r.h.s. (result in X) and compute
    !    matrix vector pruduct (result in R2: R2=A*X)
    !    workspace tranferred to prec_matv: length(R)-2*n
    !
    CALL dagmg_prec_matv(l,R,X,R(1,2),flop,R(1,3),.TRUE.)
    !
    dum(1) = DDOT(N,X,IONE,R(1,2),IONE)
    dum(2) = DDOT(N,X,IONE,R,IONE)
    BNORM = dum(3)
    rho1=dum(1)
    alpha1=dum(2)
    !
    bet0=alpha1/rho1
    !
    !
    !  Compute residual corresponding to solution bet0 * X
    CALL DAXPY(N, -bet0, R(1,2), IONE, R, IONE)
    !  ... and the corresponding square norm
    RESID = DNRM2(N,R,IONE)**2
    IF (RESID <= resi*resi*BNORM) THEN
       ! criterion satified, exit with solution bet0 * X
       ! Galerkin projection: the solution computed in this way
       ! is such that the corresponding residual
       !      R - bet0*A*X  = R - bet0*R2
       ! is orthogonal to X
       ! Indeed: bet0= ( X , R ) / ( X , R2 )
       CALL DSCAL( N, bet0, X, IONE )
       flop=flop+dble(11*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    !
    ITER = 2
    !
    !  Apply preconditioner to this new residual (result in R3) and compute
    !    matrix vector pruduct (result in R4: R4=A*R3)
    !  R3 is thus the second direction vector
    !    workspace tranferred to prec_matv: length(R)-4*n
    !
    CALL dagmg_prec_matv(l,R,R(1,3),R(1,4),flop,R(1,5),.TRUE.)
    !
    !
    dum(1) = DDOT(N,R(1,3),IONE,R(1,2),IONE)
    dum(2) = DDOT(N,R(1,3),IONE,R,IONE)
    dum(3) = DDOT(N,R(1,3),IONE,R(1,4),IONE)
    IF (.NOT.spd) dum(4) = DDOT(N,X,IONE,R(1,4),IONE)
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    IF (spd) THEN
       ! in the SPD case:
       ! gamm1 := ( X , R4 ) = ( X, A*R3 ) = ( A*X , R3 ) = ( R2 , R3)
       gamm1 = gamm0
    ELSE
       gamm1 = dum(4)
    END IF
    !
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    !
    CALL DSCAL(N, bet2, X, IONE)
    CALL DAXPY(N, zeta, R(1,3), IONE, X, IONE)
    ! Galerkin projection: the solution computed in this way
    ! is such that the corresponding residual
    !      R - bet2*A*X - zeta*A*R3 = R - bet2*R2 - zeta*R4
    ! is orthogonal to X and R3
    ! Indeed:
    !       ( X , R ) - bet2* ( X , R2 ) - zeta* ( X , R4 )
    !          =  alpha1 - bet2*rho1 - zeta*gamm1
    !          =  alpha1 - (alpha1-alpha2*gamm1/bet1) - (alpha2/bet1)*gamm1
    !          =  0
    !       ( R3 , R ) - bet2* ( R3 , R2 ) - zeta* ( R3 , R4 )
    !           [!: alpha2 is the dot product with update residual ]
    !          = ( R3 , R-bet0*R2) + bet0* ( R3 , R2)
    !                              - bet2* ( R3 , R2 ) - zeta* ( R3 , R4 )
    !          =  alpha2 + bet0*gamm0 - bet2*gamm0 - zeta*rho2
    !          =  alpha2 + alpha1*gamm0/rho1
    !                    - (alpha1-alpha2*gamm1/bet1)*(gamm0/rho1)
    !                    - (alpha2/bet1)*(bet1+gamm0*gamm1/rho1)
    !          =  0
    !
    !
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(19*N)
    IF (.NOT.spd) flop=flop+dble(2*N)
    !
    RETURN
    !
  END SUBROUTINE dagmg_Galerkin_inner
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_prolag(n, nc, V, B, ind, flop)
!
!   Prolongate vetor B based on ind:
!      P(i,j)=1 iff ind(i)=j (i=1,n ; j=1,nc)
!   Add the result to V
!
! flop : floating point ops counter
!
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    REAL(kind(0.0d0)) :: V(n), B(nc)
    REAL(kind(0.0d0)) :: flop
    !
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) THEN
          V(i) = V(i)+B(k)
       ENDIF
    ENDDO
    flop = flop + dble(n)
    RETURN
  END SUBROUTINE dagmg_prolag
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_restag(n, nc, V, B, ind, flop)
!
!   Restrict vetor V based on ind:
!      P(i,j)=1 iff ind(i)=j (i=1,n ; j=1,nc)
!   Result in B
!
! flop : floating point ops counter
!
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    REAL(kind(0.0d0)) :: V(n), B(nc)
    REAL(kind(0.0d0)) :: flop
    !
    B(1:nc)=0.0d0
    !
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) B(k) = B(k) + V(i)
    ENDDO
    flop = flop + dble(n)
    RETURN
  END SUBROUTINE dagmg_restag
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of ROUTINES FOR ITERATIVE SOLUTION
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! ROUTINES FOR SETUP
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_setupL1(n,a,ja,ia,listrank,ifl)
!
! Performs setup (corase grid definition) at level 1
! Input: matrix (n,a,ja,ia,dt(1)%idiag,dt(1)%iext)
!        in csr format with partial ordering of the rows
! Output: coarse grid matrix
!        (nc,dt(2)%a,dt(2)%ja,dt(2)%ia,dt(2)%idiag)
!        in csr format with partial ordering of the rows
!        and aggregation map in dt(1)%ind
!        - except if the level is detected to be the last one:
!           then, call to coarsest level set up -
!        The original matrix is also transformed in d+l+u format
!        (see subroutine csrdlu), stored in
!        (n,dt(1)%a,dt(1)%ja,dt(1)%il,dt(1)%iu,dt(1)%iext)
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n
    INTEGER :: ja(*),ia(n+1)
    REAL(kind(0.0d0)) :: a(*)
    INTEGER :: ifl,listrank(ifl:*)
    INTEGER :: nc,ierr,i,j,k,nz
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)) :: ops,eta,dum(2)
    CHARACTER(len=13) :: prtint
    REAL (kind(0.0d0)) :: fff(1)
    !
    nn(1)=n
    nlc(1)=n
    nlc(2)=ia(n+1)-ia(1)
    !
    ngl=nlc
       !
       wcplex=1.0d0
       nlctot=nlc
       ngltot=ngl
       ngl1=ngl
       nlc1=nlc
       icum=1
       fracnz(1)=1.0d0
       allzero=.FALSE.
       nglp=0.0d0
       maxcoarset=maxcoarsesize
       maxcoarseslowt=maxcoarsesizeslow
       maxcoarset=FLOOR(maxcoarset*(ngl(1)**(1.0d0/3)))
       maxcoarseslowt=FLOOR(maxcoarseslowt*(ngl(1)**(1.0d0/3)))
       !
       IF (wfo) THEN
          IF (n > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,'()')
          WRITE(iout,918) prtint(1:12)
          IF (nlc(2) > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,919) prtint(1:12),dble(nlc(2))/dble(n)
          WRITE(iout,'()')
       END IF
    IF ( 0 == nstep  .OR. 1 == maxlev .OR. ngl(1) <= maxcoarset ) nlev=1
    nlcp=nlc
    nglp=ngl
    !
    IF (1 /= nlev) THEN
       !
       ! Coarsening by aggregation
       !
       CALL dagmg_aggregation(1,n,a,ja,ia,nc)
       !
       ! Perform setup at next level
       !
       CALL dagmg_setup(2,nc)
       !
    ELSE
       !
       ! Coarsest level set up
       !
       DEALLOCATE(dt(1)%idiag)
       memi=memi-(n+1)
       CALL dagmg_MUMPSseq(n,a,ja,ia,fff,1,ops)
       IF (wfo) THEN
          WRITE(iout,911) IRANK,ops/dble(11*nn(1)+2*nlc1(2))
       END IF
       wcplex(4)=0.0d0
       IF (ngl(1) > maxcoarset)   wcplex(4)=-1.0d0
       IF (ngl(1) > maxcoarseslowt) wcplex(4)=ngl(1)/(1000*ngl1(1)**(1.0d0/3))
       !
       ! Work done: return
       !
       RETURN
    END IF
    !
    ! Instruction here are executed aftern return from recursive call;
    !   i.e., bottom to top
    !
       ! Change the internal structure from csr with partial ordering of the
       ! rows to separate d+l+u format (see subroutine csrdlu for details)
       !
       nz=ia(n+1)-ia(1)
       IF (transint) THEN
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%il(n+1),dt(1)%iu(n+1))
          memi=memi+nz+n+2
          memr=memr+nz
       ELSE
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%iu(n+1))
          dt(1)%il => dt(1)%idiag
          memi=memi+nz+1
          memr=memr+nz
       END IF
       memax=MAX(memax,memr+memi*rlenilen)
       !
       CALL dagmg_csrdlu(n,a,ja,ia,dt(1)%idiag              &
            ,ap,jap,dt(1)%il,dt(1)%iu,transint )
       IF (transint) THEN
          DEALLOCATE(dt(1)%idiag)
          memi=memi-(n+1)
       ELSE
          NULLIFY(dt(1)%idiag)         ! il points to same memory
       END IF
       dt(1)%a  => ap
       dt(1)%ja => jap
       NULLIFY(ap,jap)
       !
       ! Computes some internal parameter
       ! (check if K-cycle is allowed and memory needed for iterative solve)
       !
       innermax(nlev)=0
       innermax(1)=1
       eta=xsi/((1-xsi)*(cplxmax-1))
       icum=1
       DO i=2,nlev-1
          innermax(i)=min(2,floor(xsi**(i-1)/(eta*fracnz(i)*icum)))
          IF (nlev-i.LE.nlvcyc .AND. i.GT.2) innermax(i)=1
          icum=icum*innermax(i)
          wcplex(2)=wcplex(2)+icum*fracnz(i)
          wcplex(3)=wcplex(3)+(2**(i-1))*fracnz(i)
       END DO
       wcplex(2)=wcplex(2)+(2**(nlev-1))*fracnz(nlev)
       wcplex(3)=wcplex(3)+(2**(nlev-1))*fracnz(nlev)
       IF (nsmooth*smoothtype > 1)  THEN
          nwrkcum=2*nn(nlev-1)
       ELSE
          nwrkcum=nn(nlev)
       END IF
       nwrkcum=max(nwrkcum,1)
       DO i=nlev-2,1,-1
          nwrkcum=nwrkcum+3*nn(i+1)
          IF (innermax(i+1) > 1) nwrkcum=nwrkcum+2*nn(i+1)
          IF (nsmooth*smoothtype > 1)  nwrkcum=max(2*nn(i),nwrkcum)
       END DO
       IF (wfo) THEN
          WRITE(iout,'()')
          WRITE(iout,954) nlctot(1)/dble(nlc1(1))
          WRITE(iout,955) nlctot(2)/dble(nlc1(2))
       END IF
       IF (wff) THEN
          WRITE(iout,956) wcplex(3)
          WRITE(iout,957) wcplex(2)
          WRITE(iout,'()')
       END IF
    RETURN
911 FORMAT(i3,'*','        Exact factorization:',f12.3,            &
         ' equiv. CG iterations')
918 FORMAT('****','        Number of unknowns:', A12)
919 FORMAT('****','                 Nonzeros :', A12,                &
         ' (per row:',f7.2,')')
920 FORMAT('****','        Number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
921 FORMAT('****','                   Nonzeros:',A12,              &
         ' (per row:',f4.1,    &
         '; red. ratio:',f5.2,')')
954 FORMAT('****','                  Grid complexity:',f9.2)
955 FORMAT('****','              Operator complexity:',f9.2)
956 FORMAT('****','  Theoretical Weighted complexity:',f9.2, &
                ' (K-cycle at each level)' )
957 FORMAT('****','    Effective Weighted complexity:',f9.2, &
                ' (V-cycle enforced where needed)' )
  END SUBROUTINE dagmg_setupL1
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE dagmg_setup(l,n,listrank)
!
! Performs setup (corase grid definition) at level l > 1
! Input: matrix (n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag,dt(l)%iext)
!        in csr format with partial ordering of the rows
! Output: coarse grid matrix
!        (nc,dt(l+1)%a,dt(l+1)%ja,dt(l+1)%ia,dt(l+1)%idiag,dt(l+1)%iext)
!        in csr format with partial ordering of the rows
!        and aggregation map in dt(l)%ind
!        - except if the level is detected to be the last one:
!           then, call to coarsest level set up -
!        The original matrix is also transformed in d+l+u format
!        (see subroutine csrdlu), stored in
!        (n,dt(l)%a,dt(l)%ja,dt(l)%il,dt(l)%iu,dt(l)%iext)
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n
    INTEGER, OPTIONAL :: listrank(n+1:*)
    INTEGER :: nc,ierr,i,j,k,nz
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    LOGICAL, SAVE :: slowcoarse
    REAL(kind(0.0d0)) :: ops,fw,eta,dum(2)
    CHARACTER(len=13) :: prtint
    REAL (kind(0.0d0)) :: fff(1)
    !
    nn(l)=n
    nlc(1)=n
    IF (n > 0) THEN
       nlc(2)=dt(l)%ia(n+1)-dt(l)%ia(1)
    ELSE
       nlc(2)=0
    END IF
    !
    ngl=nlc
       !
       nlctot=nlctot+nlc
       ngltot=ngltot+ngl
       IF (wfo) THEN
          WRITE(iout,'()')
          WRITE(iout,914) IRANK,l
          IF (n > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,920) prtint(1:12),dble(nlcp(1))/dble(n)
          IF (nlc(2) > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,921) prtint(1:12),dble(nlc(2))/dble(n),        &
               dble(nlcp(2))/dble(nlc(2))
       END IF
       fracnz(l)=ngl(2)/ngl1(2)
    IF (l==2) THEN
       wcplex(1)=ngl1(2)/ngl(2)
       slowcoarse=.FALSE.
    END IF
       !
    slcoarse = nglp(1) < 2*ngl(1) .AND. nglp(2) < 2*ngl(2)
    IF( l == nstep+1  .OR. l == maxlev                   &
         .OR. (  ngl(1) <= maxcoarset)                   &
         .OR. ( slcoarse .AND. ngl(1) <= maxcoarseslowt )  &
!         .OR. ( zerors .AND. ngl(1) <= maxcoarseslowt )  &
         .OR. ( slowcoarse .AND. slcoarse )              &
         .OR. nglp(1) ==ngl(1) )                       THEN
       nlev=l
    END IF
    slowcoarse=slcoarse
    nlcp=nlc
    nglp=ngl
    IF (n == 0) THEN
       nc=0
       allzero=.TRUE.
       !
       ! Work done: return
       !
       RETURN
    END IF
    IF (l /= nlev) THEN
       !
       ! Coarsening by aggregation
       !
       CALL dagmg_aggregation(l,n,dt(l)%a,dt(l)%ja,dt(l)%ia,nc)
       !
       ! Perform setup at next level
       !
       CALL dagmg_setup(l+1,nc)
       !
    ELSE
       !
       ! Coarsest level set up
       !
       DEALLOCATE(dt(l)%idiag)
       memi=memi-(n+1)
       CALL dagmg_MUMPSseq(n,dt(l)%a,dt(l)%ja,dt(l)%ia,fff,1,ops)
       !
       ! Coarse  solver uses its own internal memory: a,ja,ia
       ! more needed
       memi=memi-size(dt(l)%ja)-n-1
       memr=memr-size(dt(l)%a)
       DEALLOCATE(dt(l)%a,dt(l)%ja,dt(l)%ia)
       !
       IF (wfo) THEN
          WRITE(iout,911) IRANK,ops/dble(11*nn(1)+2*nlc1(2))
       END IF
       wcplex(4)=0.0d0
       IF (ngl(1) > maxcoarset)   wcplex(4)=-1.0d0
       IF (ngl(1) > maxcoarseslowt) wcplex(4)=ngl(1)/(1000*ngl1(1)**(1.0d0/3))
       !
       ! Work done: return
       !
       RETURN
    END IF
    !
    ! Instruction here are executed aftern return from recursive call;
    !   i.e., bottom to top
    !
       ! Change the internal structure from csr with partial ordering of the
       ! rows to separate d+l+u format (see subroutine csrdlu for details)
       !
       nz=dt(l)%ia(n+1)-dt(l)%ia(1)
       IF (transint) THEN
          ALLOCATE(ap(nz),jap(nz-n),dt(l)%il(n+1),dt(l)%iu(n+1))
          memi=memi+nz+n+2
          memr=memr+nz
       ELSE
          ALLOCATE(ap(nz),jap(nz-n))
          dt(l)%iu => dt(l)%ia
          dt(l)%il => dt(l)%idiag
          memi=memi+nz-n
          memr=memr+nz
       END IF
       memax=MAX(memax,memr+memi*rlenilen)
       !
       CALL dagmg_csrdlu(n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag &
            ,ap,jap,dt(l)%il,dt(l)%iu,transint )
       memi=memi-size(dt(l)%ja)
       memr=memr-size(dt(l)%a)
       DEALLOCATE(dt(l)%a,dt(l)%ja)
       IF (transint) THEN
          DEALLOCATE(dt(l)%idiag,dt(l)%ia)
          memi=memi-2*(n+1)
       ELSE
          NULLIFY(dt(l)%idiag,dt(l)%ia)     ! il,iu point to same memory
       END IF
       dt(l)%a  => ap
       dt(l)%ja => jap
       NULLIFY(ap,jap)
    !
    RETURN
911 FORMAT(i3,'*','        Exact factorization:',f12.3,            &
         ' equiv. CG iterations')
914 FORMAT(i3,'*','                      Level:',I12)
918 FORMAT('****','        Number of unknowns:', A12)
919 FORMAT('****','                 Nonzeros :', A12,                &
         ' (per row:',f7.2,')')
920 FORMAT('****','        Number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
921 FORMAT('****','                   Nonzeros:',A12,              &
         ' (per row:',f4.1,    &
         '; red. ratio:',f5.2,')')
954 FORMAT('****',' relative complexity (grid):',f9.2)
955 FORMAT('****',' relative complexity (nnzs):',f9.2)
  END SUBROUTINE dagmg_setup
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_smoothsetup
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,i
    !
    ! Smoother set up
    !
    IF (smoothtype == 1) THEN
       DO l=1,nlev-1
          !
          ! Gauss-Seidel smoother: store inverse of diagonal elements
          ! in dt(l)%a(1:n)
          DO i=1,nn(l)
             dt(l)%a(i)=1.0d0/dt(l)%a(i)
          END DO
       END DO
    END IF
    RETURN
  END SUBROUTINE dagmg_smoothsetup
!------------------------------------------------------------------
  SUBROUTINE dagmg_aggregation(l,n,a,ja,ia,nc,listrank)
!
!  Perform aggregation of level l matrix with partial
!  ordering of the rows specified in (n,a,ja,ia,idiag,iext)
!  (iext significant only with the parallel version)
!
!  Output:
!     nc : number of aggregates (= number of coarse unknowns)
!     dt(l)%ind(1:n): ind(i) is the index of the aggregates to which
!                     i belongs; ind(i)=0 iff i has been kept outside
!                     aggregation (see checkdd).
!     P: associated prologation matrix;
!        P(i,j)=1 iff ind(i)=j and P(i,j)=0 otherwise (i=1,n ; j=1,nc).
!
!     Corresponding coarse grid matrix stored with partially ordered rows
!       in dt(l+1)%a, dt(l+1)%ja, dt(l+1)%ia, dt(l+1)%idiag, dt(l+1)%iext
!                                   (both these latter with length nc+1).
!
!  note: storage not optimal: dt(l+1)%a, dt(l+1)%ja greater than needed
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,nc
    INTEGER :: ja(*),ia(n+1)
    INTEGER, OPTIONAL :: listrank(*)
    REAL (kind(0.0d0)) :: a(*), dum
    INTEGER :: ier,i,j,k,maxdg,np,kpass,nzc,m1,ndd,nzp,isize,nddp,npass1,nz
    LOGICAL :: skipass
    !
    INTEGER, POINTER, DIMENSION(:) :: jan,ian,idiagn,iextn,ind2,lcg,lcgn
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: an
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sinn
    !
    INTEGER, POINTER, DIMENSION(:) :: jap,iap,idiagp,iextp,lcg1
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sip
    !
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ldd,iw,iperm,riperm
    REAL(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: si1,w
    REAL(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: wc
    !
    !
    IF (l .EQ. 1) THEN
       IF (wfo) THEN
          WRITE (iout,901) IRANK
       END IF
       IF (wff) THEN
          IF (spd) THEN
             WRITE (iout,906)
          ELSE IF (  transint) THEN
             WRITE (iout,908)
          END IF
          IF (.not.spd) then
             WRITE (iout,902) 'Jacobi',kaptg_dampJac,checkddJ
          ELSE
             WRITE (iout,902) 'BlockD',kaptg_blocdia,checkddB
          END IF
          IF (checkdd < 0) THEN
             WRITE (iout,904)
          END IF
          WRITE(iout,903) npass,targetcoarsefac
          WRITE (iout,905) trspos
       END IF
    END IF
    !
    IF (l .EQ. 1) THEN
       ALLOCATE(si1(n),ind2(n),iperm(n),riperm(n))
       memi=memi+4*n+1
       memr=memr+n
       CALL dagmg_setCMK(n,ja,ia,dt(l)%idiag,riperm,iperm)
    ELSE
       ALLOCATE(si1(n),ind2(n),iperm(n))
       memi=memi+2*n
       memr=memr+n
       iperm(1:n)=1
    END IF
    memax=MAX(memax,memr+memi*rlenilen)
    !
    call dagmg_prepareagg(n,a,ja,ia,dt(l)%idiag,ind2,iperm,si1,ndd,l )
    !
    IF (ndd .EQ. n) THEN
       nc=0
       nzc=0
       DEALLOCATE(si1,iperm,ind2)
       memi=memi-2*n
       memr=memr-n
       IF (l .EQ. 1) THEN
          DEALLOCATE(riperm)
          memi=memi-n
       END IF
       GOTO 999
    END IF
    !
    ALLOCATE(ldd(ndd),lcg(2*(n-ndd)))
    memi=memi+2*n-ndd
    memax=MAX(memax,memr+memi*rlenilen)
    !
    IF (dble(n) .GT. targetcoarsefac*(n-ndd)) THEN
       !  skip first pairwise pass if the coarsening seems sufficiently
       !  fast based on ndd only
       skipass=.TRUE.
       !  but add one more pass to compensate in case of need
       npass1=npass+1
    ELSE
       skipass=.FALSE.
       npass1=npass
    END IF
    !
    !  perform initial pairwise aggregation; output in ind2 and lcg
    !
    IF (l > 1) THEN
       IF (spd) THEN
          CALL dagmg_findpairs_SI(n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       ELSE
          CALL dagmg_findpairs_GI(n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       END IF
       DEALLOCATE(iperm)
       memi=memi-n
    ELSE
       IF (spd) THEN
          CALL dagmg_findpairs_SI1(n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       ELSE
          CALL dagmg_findpairs_GI1(n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       END IF
       DEALLOCATE(iperm,riperm)
       memi=memi-2*n
    END IF
    nz=ia(n+1)-1
    !
    ! Form aggregated matrix in (an,jan,ian)
    ! the new matrix has at most the following number of nonzero:
    !  nz (the previous number) - ndd - 2*(n-ndd-nc)
    !         indeed: from rows listed in ldd, at least ndd
    !         diagonal entries are removed;
    !         (n-ndd-nc) is the number of pairs formed,
    !         and each of them condensate at least 3 entries in a
    !         single diagonal entry
    IF (npass1.GT.1) THEN
       isize=nc
    ELSE
       isize=1
    END IF
    ALLOCATE(an(nz-2*(n-nc)+ndd),jan(nz-2*(n-nc)+ndd)             &
         ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
    memi=memi+nz-2*(n-nc)+ndd+2*(nc+1)+2*nc+isize
    memr=memr+nz-2*(n-nc)+ndd+nc
    memax=MAX(memax,memr+memi*rlenilen)
    CALL dagmg_setcg(n,a,ja,ia,dt(l)%idiag,si1,ind2,lcg,nc,an,jan,ian      &
         ,idiagn,sinn,npass1.GT.1,maxdg,iw,wc )
    DEALLOCATE(wc,iw)
    memi=memi-2*nc
    memr=memr-nc
    nzc=ian(nc+1)-1
    IF (dble(nz).GT.targetcoarsefac*nzc .OR. npass1.LE.1) THEN
       DEALLOCATE(si1,sinn,lcg)
       memi=memi-2*(n-ndd)
       memr=memr-n-isize
       dt(l)%ind  => ind2
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(ind2,an,jan,ian,idiagn)
       GOTO 999
    END IF
    !
    DEALLOCATE(ind2)
    memi=memi-n
    !
    lcg1 => lcg
    NULLIFY(lcg)
    m1=1
    !
    !
    ! Allocated pointers at this stage: an, jan, ian, idiagn, sinn, lcg1
    !                                             (lcg accessed via lcg1)
    ! Allocated arrays: si1, ldd
    !
    DO kpass=2,npass1
       m1=2*m1
       np  = nc
       nzp = nzc
       ap     => an
       jap    => jan
       iap    => ian
       idiagp => idiagn
       sip    => sinn
       NULLIFY(an,jan,ian,idiagn,sinn)
       !
       ALLOCATE(lcg(2*np),ind2(np),w(maxdg),iw(maxdg))
       memi=memi+maxdg+3*np
       memr=memr+maxdg
       memax=MAX(memax,memr+memi*rlenilen)
       !
       !   perform further pairwise aggregation; output in ind2 and lcg
       !
       ind2(1:np)=-1
       IF (spd) THEN
          CALL dagmg_findpairs_SF(np,ap,jap,iap,idiagp,sip  &
               ,ind2,lcg,nc                                  &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       ELSE
          CALL dagmg_findpairs_GF(np,ap,jap,iap,idiagp,sip     &
               ,ind2,lcg,nc                                     &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       END IF
       DEALLOCATE(w,iw)
       memi=memi-maxdg
       memr=memr-maxdg
       IF (kpass.NE.npass1) THEN
          isize=nc
       ELSE
          isize=1
          DEALLOCATE(si1)
          memr=memr-n
       END IF
       !
       ! Form aggregated matrix in (an,jan,ian,idiagn,iextn)
       ! the new matrix has at most the following number of nonzero:
       !  nzp (the previous number) - 2*(np-nc)
       !         indeed: np-nc is the number of pairs formed,
       !         and each of them condensate at least 3 entries in a
       !         single diagonal entry
       !
       !  below we reallocate an,ja,ian,idiagn; the memory they were
       !        pointing to can be accessed via ap,jap,iap,idiagp
       ALLOCATE(an(nzp-2*(np-nc)),jan(nzp-2*(np-nc))                 &
            ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
       memi=memi+nzp-2*(np-nc)+2*(nc+1)+2*nc+isize
       memr=memr+nzp-2*(np-nc)+nc
       memax=MAX(memax,memr+memi*rlenilen)
       !
       ! Allocated pointers at this stage: an,jan,ian,idiagn,sinn,(iextn)
       !                                   ap,jap,iap,idiagp,sip,(iextp)
       !                                   lcg1,lcg,ind2
       ! Allocated arrays: si1,ldd,wc,iw
       !
       CALL dagmg_setcg(np,ap,jap,iap,idiagp,sip,ind2,lcg,nc,an     &
            ,jan,ian,idiagn,sinn,kpass.NE.npass1,maxdg,iw,wc )
       memi=memi-SIZE(jap)-2*(np+1)-np-2*nc
       memr=memr-SIZE(ap)-SIZE(sip)-nc
       DEALLOCATE(ap,jap,iap,idiagp,sip,ind2,wc,iw)
       !
       ! Allocated pointers at this stage: an,jan,ian,idiagn,(iextn)
       !                                   sinn,lcg1,lcg
       ! Allocated arrays: si1,ldd
       !
       ! Form new list of aggregates from previou one and the pairwise step
       !   just performed
       ALLOCATE(lcgn(2*m1*nc))
       memi=memi+2*m1*nc
       memax=MAX(memax,memr+memi*rlenilen)
       CALL dagmg_lcgmix(nc,m1,lcg1,lcg,lcgn)
       memi=memi-SIZE(lcg)-SIZE(lcg1)
       DEALLOCATE(lcg,lcg1)
       lcg1 => lcgn
       NULLIFY(lcgn)
       !
       ! Allocated pointers at this stage: an,jan,ian,idiagn,sinn,lcg1,(iextn)
       ! Allocated arrays: si1,ldd
       !  -- same as before entering the do loop --
       nzc=ian(nc+1)-1
       IF ( kpass.NE.npass1 .AND. dble(nz).GT.targetcoarsefac*nzc ) THEN
          DEALLOCATE(si1)
          memr=memr-n
          EXIT
       END IF
    END DO
    !
    memr=memr-SIZE(sinn)
    DEALLOCATE(sinn)
    !
    ALLOCATE(dt(l)%ind(n))
    memi=memi+n
    memax=MAX(memax,memr+memi*rlenilen)
    CALL dagmg_setind(nc,ndd,ldd,lcg1,2*m1,dt(l)%ind)
    memi=memi-ndd-SIZE(lcg1)
    DEALLOCATE(lcg1,ldd)
    !
    ! Allocated pointers at this stage: an,jan,ian,idiagn,(iextn)
    ! Allocated arrays: -
    !
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(an,jan,ian,idiagn)
    !
999 CONTINUE
    !
    RETURN
901 FORMAT(i3,'*SETUP: Coarsening by multiple pairwise aggregation')
902 FORMAT('****       Quality threshold (',A6,'):',f6.2, &
         ' ;  Strong diag. dom. trs:',f5.2)
903 FORMAT('****         Maximal number of passes:',i3,     &
         '  ; Target coarsening factor:',f5.2)
904 FORMAT('****           Diag. dom. checked w.r.t. sum of offdiag', &
          ' (no absolute vaues)')
905 FORMAT('****',22x,'Threshold for rows with large pos. offdiag.:',f5.2)
906 FORMAT('****  Rmk: Setup performed assuming the matrix symmetric')
908 FORMAT('****  Rmk: Setup performed for the transpose of the input matrix')
  END SUBROUTINE dagmg_aggregation
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_setCMK(n,ja,ia,idiag,riperm,iperm,iext)
    !
    !     compute CMK permutation
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),riperm(*),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    LOGICAL :: exc
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,mindg,kdim
    !
    !
    ! First find starting node of minimal degree
    ! while already numbering nodes whith degree one
    ! (i.e., only the diagonal element is present in the row);
    ! store the opposite of the degree in iperm
    mindg=n+1
    jj=1
    i2=1
    DO i = 1,n
       dg=ia(i+1)-ia(i)
       IF (dg .GT. 1) THEN
          iperm(i)=-dg
          IF (dg.LT.mindg) THEN
             mindg=dg
             jj=i
          END IF
       ELSE
          riperm(i2)=i
          iperm(i)=i2
          ! one more node is numbered: increase i2
          i2=i2+1
       END IF
    ENDDO
    !
    ijs=0
    i1=i2
15  CONTINUE
    !
    ! Start with next node of minimum degree: jj
    riperm(i2)=jj
    iperm(jj)=i2
    !
    ! in what follows: i2 is the number of nodes already numbered
    !    and i1-1 the number of nodes whose neighbors have been processed
    DO WHILE (i1.LE.i2 .AND. i2.LT.n)
       !
       ! Number neighbors of riperm(i1)
       !
       i=riperm(i1)
       ijs1=i2+1
       !
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       DO kk = j1,jd-1
          j=ja(kk)
          IF (iperm(j) .LT. 0) THEN
             ! one more node is numbered: increase i2
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       DO kk = jd+1,j2
          j=ja(kk)
          IF (iperm(j) .LT. 0) THEN
             ! one more node is numbered: increase i2
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       !
       ijs2=i2
       ! Sort just found nodes such that value stored in iperm
       !  will be in descending order; since this is the oposite of
       !  the degree, degrees will be in ascending order
       exc=.TRUE. .AND. ijs2.GT.ijs1
       DO WHILE(exc)
          exc=.FALSE.
          DO kk=ijs1+1,ijs2
             IF( iperm(riperm(kk)) .GT. iperm(riperm(kk-1)) )THEN
                j=riperm(kk)
                riperm(kk)=riperm(kk-1)
                riperm(kk-1)=j
                exc=.TRUE.
             END IF
          END DO
       END DO
       DO kk=ijs1,ijs2
          iperm(riperm(kk))=kk
       END DO
       !
       ! Done with riperm(i1): proceed with next node
       i1=i1+1
    END DO
    IF (i2 .LT. n) THEN
       !
       ! Not all nodes numbered during first pass
       !   (may happen if the matrix is reducible)
       ! Find a new starting node with minimal degree and proceed
       jj=0
       DO WHILE (jj .EQ. 0)
          ijs=ijs+1
          IF (ijs .GT. n) THEN
             mindg=mindg+1
             ijs=1
          END IF
          ijs1=ijs
          IF (iperm(ijs1).LT.0 .AND. ia(ijs1+1)-ia(ijs1).EQ.mindg) &
               jj=ijs1
       END DO
       i2=i2+1
       ! Repeat previous step
       GOTO 15
    END IF
    !
    RETURN
  END SUBROUTINE dagmg_setCMK
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_prepareagg(n,a,ja,ia,idiag,ind2,iperm,si,ndd,l,iext)
    !
    !  Performs some setup before aggregation:
    !     detects nodes to be kept outside because of strong
    !     diagonal dominance, as well as nodes to be transferred
    !     as is to the coarse grid because of strong positive
    !     coupling in the row or column
    !
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind2(n),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER :: ndd, l
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) , TARGET :: si(n)
    REAL(kind(0.0d0)) :: checkddl,oda,odm,ods,vald
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,kdim,nnegrcs
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: odmax,odabs,osi
    !
    IF (.NOT.spd) THEN
       checkddl=checkddJ
    ELSE
       checkddl=checkddB
    END IF
    !
    IF (.NOT.spd) THEN
       !
       ! Matrix is not symmetric: gather for each node
       !   information about elements in column:
       !   sum, maximum and, if needed, sum of abs value
       !
       IF (checkdd > 0) THEN
          ALLOCATE(odmax(n),odabs(n))
          memi=memi+2*n
          odabs(1:n)=0.0d0
       ELSE
          ALLOCATE(odmax(n))
          memi=memi+n
       END IF
       osi => si(1:n)
       si(1:n)=0.0d0
       odmax(1:n)=0.0d0
       memax=MAX(memax,memr+memi*rlenilen)
       DO i=n,1,-1
          j =ia(i)
          jd=idiag(i)
          jj=ia(i+1)-1
          DO k = j,jd-1
             jk=ja(k)
             osi(jk)=osi(jk)+a(k)
             odmax(jk)=max(odmax(jk),a(k))
             IF (checkdd > 0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
          DO k = jd+1,jj
             jk=ja(k)
             osi(jk)=osi(jk)+a(k)
             odmax(jk)=max(odmax(jk),a(k))
             IF (checkdd > 0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
       ENDDO
    END IF
    !
    ! Compute the mean row and column sum,
    !  and check if the node will participate to aggregation
    !
    ndd=0
    nnegrcs=0
    !
    DO i=1,n
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       vald = a(jd)
       odm=0.0d0
       oda=0.0d0
       ods=0.0d0
       DO kk = j1,jd-1
          ods=ods+a(kk)
          odm=max(odm,a(kk))
          IF (checkdd > 0) oda=oda+abs(a(kk))
       ENDDO
       DO kk = jd+1,j2
          ods=ods+a(kk)
          odm=max(odm,a(kk))
          IF (checkdd > 0) oda=oda+abs(a(kk))
       ENDDO
       !
       ! Update with info from columns if needed
       IF (.NOT.spd) THEN
          ods=(osi(i)+ods)/2
          odm=max(odm,odmax(i))
          IF (checkdd > 0) oda=(oda+odabs(i))/2
       END IF
       !
       ! Count the number of nodes with negative mean row and column sum
       IF ((vald+ods) .LT. -repsmach*ABS(vald)) nnegrcs=nnegrcs+1
       !
       ! Now check the category of the node and fill si(i)
       !
       si(i)=-ods
       IF ( (checkdd.GT.0 .AND. vald.GT.checkddl*oda)     &
            .OR. (checkdd.LT.0 .AND. vald.GT.checkddl*abs(ods)) ) THEN
          !
          ! Node to be kept outside the aggregation because of strong diagonal
          !   dominance: set ind2(i) to 0
          ind2(i)=0
          ndd=ndd+1
       ELSE
          ! Normal node: set ind2(i) to negative
          ind2(i)=-1
          !
          ! iperm(i) set to zero for nodes that are to be transeferred as is
          ! to the coarse grid
          IF (odm .GT. trspos*vald) iperm(i)=0
       ENDIF
    END DO
    !
    ! set zerors according to nnegrcs and fracnegrcsum
    !
    zerors=.FALSE.
    IF (nnegrcs .GT. fracnegrcsum*n) THEN
       zerors=.TRUE.
       ndd=0
       ind2(1:n)=-1
    END IF
    !
    ! Release local memory if needed
    IF (.NOT.spd) THEN
       IF (checkdd > 0) THEN
          DEALLOCATE(odmax,odabs)
          memi=memi-2*n
       ELSE
          DEALLOCATE(odmax)
          memi=memi-n
       END IF
    END IF
    RETURN
  END SUBROUTINE dagmg_prepareagg
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_GF(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
! Performs pairwise aggregation according to Algorithms 4.2 and 4.3 in [2,3],
! with some heuristic enhancements to cover non diagonally  dominant matrices.
!
! Version: further pairwise aggregation for general matrices [3].
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    REAL(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
!
!  INPUT: n,ja(*),ia(n+1),idiag(n): aggregated matrix from previous pass
!           (CSR format with partial ordering of each row: lower
!            triangular part, next diag., and eventually upper triagular part;
!            idiag point to the diagonal element)
!
!         si(n) : vector s or \tile{s} in Algorithm 4.2 or 4.3 of [2,3]
!
!         ind(n): on input, entries in ind should be negative
!
!         ja1(*),ia1(*),idiag1(*): unaggregated matrix, in which the
!                threshold condition has to be checked
!         si1(*): si in that matrix
!
!         m1    : maximal number of nodes in aggregates from previous passes
!
!         lcg1(m1,n): list of aggregates from previous passes
!           lcg1(2,i) = 0 means that node i has to be transferred unaggregated
!                        to the coarse grid, and lcg(2,j) set to 0
!                        for the corresponding corase grid node j
!
!         zerors: if true, diagonal entries are not taken into account;
!                 instead, the matrix is treated as it would have mean
!                 column and row sum equal to zero
!
!  OUTPUT: nc   : number of pairs formed
!
!          lcg(2,nc): list of aggregates;
!                     lcg(2,i)=0 : lcg(1,i) is a singleton that was
!                                  forced to be transferred unaggregated
!                     lcg(2,i)=-1: lcg(1,i) is a singleton because no valid
!                                  pair has been found
!                     lcg(2,i)>0 : (lcg(1,i),lcg(2,i)) form the ith pair
!
!          ind(n): ind(i)=j means that i belongs to aggregate j
!
!  rtent, jtent are working arrays of size >= the maximal degree of a node in
!               input matrix (a,ja,ia)
!
!----------------
! Local variables
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       ijs=ijs+1
       !
       ! Node isel has been selected
       !
       ! Check if isel has already been processed
       IF (ind(isel) .GE. 0) CYCLE
       !
       ! A new aggregate is formed that contains isel
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       ! Check if isel has to be transferred unaggregated
       IF (lcg1(2,isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       !
       !  Try to form a pair with isel: follow list of neighbors
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          !  CYCLE means: reject this neighbor, proceed with next one
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          !  check if j is available to form a pair
          IF(lcg1(2,j).EQ.0 .OR. ind(j).GE.0) CYCLE
          !   search for the corresponding entry in transposed matrix
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !   CYCLE instructions below: pair rejected because A_G is not
          !   nonnegative definite
          !
          !  Heuristic if sigj is negative: set Sigma_G(j)=abs(sigj) (j=1,2)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          !   main threshold test
          IF (valp > kaptg) CYCLE
          !
          !    A_G is nonneagtive definite and passed the corresponding
          !    "quality" threshold test: (isel,j) is possibly
          !    an acceptable pair; check if it is the best one
          !      and update the list of possible second choices in rtent, jtent
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair .EQ. 0) GOTO 25
20     CONTINUE
       !       Perform check in unagregated matrix, considering possible pairs
       !       in right order
       CALL dagmg_checktentagg_GF
       IF (.NOT.acc) THEN
          !     Tentative pair has been rejected, check for the next one, if any
          ipair = 0
          IF (ntentleft .GT.0) THEN
             i=1
             j=1
             DO WHILE (i .LE. ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
       !
25     CONTINUE
       !
       IF (ipair .EQ. 0) THEN
       ! no valid pair found: isel left alone in aggregate nc
          lcg(2,nc) = -1
       ELSE
       ! pair nc is formed with isel and ipair
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
    RETURN
  CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_checktentagg_GF
!
!  Check the quality in the original matrix of an aggregate formed
!   by grouping two aggregates from previous pass(es)
!
!  INPUT: see dagmg_findpairs_GF
!
!  OUTPUT: acc(logical): acc=.true. if the new aggregate can be accepted and
!                        acc=.false. otherwise
!
!
    INTEGER, PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, f1, f2
    INTEGER :: j,jj,k,l,m,info, setdim1, setdim, l2, k2
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
!
! Find indices of submatrix to be extracted and store them in set(1:setdim)
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel) .LT. 0) THEN
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
!
! Sort indices in set(1:setdim) in increasing order
!   (makes susequent processing faster)
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l)<set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
!
! Extract submatrix of (the symmetric part of) input matrix and store it in W
! Store in AGe the mean row and column sum of input matrix for related indices
! (this is also the row sum of A_G)
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=a1(idiag1(jj))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)/2
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
       DO k=ia1(jj),idiag1(jj)-1
          DO l=1,j-1
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)/2
                W(j,l)=W(j,l)+alpha
                W(l,j)=W(j,l)
                EXIT
             END IF
          END DO
       END DO
    END DO
!
    DO j=1,SetDim
   !       Set sig(j) equal to minus the sum of connections that are
   !       external to the tentative aggregate; that is, add the internal
   !       offdiagonal entries to current sig(j)
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
   !       Heuristic if sig(j) is negative: set Sigma_G(j)=abs(sig(j))
   !         Hence need to correct the row sum of A_G stored in AGe
       IF (sig(j) < 0.0d0)  AGe(j)=AGe(j)+2*sig(j)
   !
   !       Store diagonal element in v
       v(j)=W(j,j)
   !                          2
   !       Now set W = A_G - --- D_G
   !                         kap
   !
   !       Then W contains the matrix that need to be nonnegative definite to
   !       accept the aggregate according to the theory in [3],
   !       up to the rank one correction that is added below
   !
   !        Offdiagonal entrie of W alraedy OK;
   !        for the diagonal ones, note that A_G(j,j)=W(j,j)-abs(sig(j))
   !        and D_G(j,j)=W(j,j), hence the rule below with umdbndmum1=1-2/kap
   !
       W(j,j)=umdbndmum1*W(j,j)-abs(sig(j))
   !       beta <- quadratic form (1^T * D_G * 1) ; [ 1 = vector of all ones ]
   !       alp is positive  if and only if AGe(j) is zero for all j;
   !              that is, if and only if A_G has 1 in its null space
       IF (j .eq. 1) THEN
          beta=v(j)
          alp=abs(AGe(j))
       ELSE
          beta=beta+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
!
! Eventually, add to W the rank one correction
!                2
!       -------------------  D_G * 1 * 1^T * D_G
!       kap*(1^T * D_G * 1)
!
    beta=dbndmum1/beta
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
!
! Now, the rule is to accept the tentative aggregate if and only if
! W(setdim,setdim) is nonnegative definite
! Two cases:
!     alp == 0 , then W is singular, hence we check the positive
!                definitenes of W(setdim-1,setdim-1)
!     alp > 0  , then W is normally not singular and we require thus that
!                it is positive definite
!
    IF (alp.LT.repsmach*beta) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
!
! Perform a Choleky factorization of W(setdim1,setdim1) and check
! for negative pivots; code expended for small dimensions (usual case)
!
    acc=.FALSE.
! acc has just be defined to FALSE; hence RETURN statement means rejection
!
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info .NE. 0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8) .LE. 0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7) .LE. 0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6) .LE. 0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5) .LE. 0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4) .LE. 0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3) .LE. 0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2) .LE. 0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1) .LE. 0.0d0) RETURN
10  CONTINUE
!
!   all test passed: accept the aggregate
    acc=.TRUE.
!
    RETURN
  END SUBROUTINE dagmg_checktentagg_GF
  END SUBROUTINE dagmg_findpairs_GF
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_SF(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
! Performs pairwise aggregation according to Algorithms 4.2 and 4.3 in [2,3],
! with some heuristic enhancements to cover non diagonally  dominant matrices.
!
! Version: further pairwise aggregation for symm. matrices [2].
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    REAL(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
!
!  INPUT: n,ja(*),ia(n+1),idiag(n): aggregated matrix from previous pass
!           (CSR format with partial ordering of each row: lower
!            triangular part, next diag., and eventually upper triagular part;
!            idiag point to the diagonal element)
!
!         si(n) : vector s or \tile{s} in Algorithm 4.2 or 4.3 of [2,3]
!
!         ind(n): on input, entries in ind should be negative
!
!         ja1(*),ia1(*),idiag1(*): unaggregated matrix, in which the
!                threshold condition has to be checked
!         si1(*): si in that matrix
!
!         m1    : maximal number of nodes in aggregates from previous passes
!
!         lcg1(m1,n): list of aggregates from previous passes
!           lcg1(2,i) = 0 means that node i has to be transferred unaggregated
!                        to the coarse grid, and lcg(2,j) set to 0
!                        for the corresponding corase grid node j
!
!         zerors: if true, diagonal entries are not taken into account;
!                 instead, the matrix is treated as it would have mean
!                 column and row sum equal to zero
!
!  OUTPUT: nc   : number of pairs formed
!
!          lcg(2,nc): list of aggregates;
!                     lcg(2,i)=0 : lcg(1,i) is a singleton that was
!                                  forced to be transferred unaggregated
!                     lcg(2,i)=-1: lcg(1,i) is a singleton because no valid
!                                  pair has been found
!                     lcg(2,i)>0 : (lcg(1,i),lcg(2,i)) form the ith pair
!
!          ind(n): ind(i)=j means that i belongs to aggregate j
!
!  rtent, jtent are working arrays of size >= the maximal degree of a node in
!               input matrix (a,ja,ia)
!
!----------------
! Local variables
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       ijs=ijs+1
       !
       ! Node isel has been selected
       !
       ! Check if isel has already been processed
       IF (ind(isel) .GE. 0) CYCLE
       !
       ! A new aggregate is formed that contains isel
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       ! Check if isel has to be transferred unaggregated
       IF (lcg1(2,isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       !
       !  Try to form a pair with isel: follow list of neighbors
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          !  CYCLE means: reject this neighbor, proceed with next one
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          !  check if j is available to form a pair
          IF(lcg1(2,j).EQ.0 .OR. ind(j).GE.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !   CYCLE instructions below: pair rejected because A_G is not
          !   nonnegative definite
          !
          !  Heuristic if sigj is negative: set Sigma_G(j)=abs(sigj) (j=1,2)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          !   main threshold test
          IF (valp > kaptg) CYCLE
          !
          !    A_G is nonneagtive definite and passed the corresponding
          !    "quality" threshold test: (isel,j) is possibly
          !    an acceptable pair; check if it is the best one
          !      and update the list of possible second choices in rtent, jtent
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair .EQ. 0) GOTO 25
20     CONTINUE
       !       Perform check in unagregated matrix, considering possible pairs
       !       in right order
       CALL dagmg_checktentagg_SF
       IF (.NOT.acc) THEN
          !     Tentative pair has been rejected, check for the next one, if any
          ipair = 0
          IF (ntentleft .GT.0) THEN
             i=1
             j=1
             DO WHILE (i .LE. ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
       !
25     CONTINUE
       !
       IF (ipair .EQ. 0) THEN
       ! no valid pair found: isel left alone in aggregate nc
          lcg(2,nc) = -1
       ELSE
       ! pair nc is formed with isel and ipair
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
    RETURN
  CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_checktentagg_SF
!
!  Check the quality in the original matrix of an aggregate formed
!   by grouping two aggregates from previous pass(es)
!
!  INPUT: see dagmg_findpairs_SF
!
!  OUTPUT: acc(logical): acc=.true. if the new aggregate can be accepted and
!                        acc=.false. otherwise
!
!
    INTEGER, PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, f1, f2
    INTEGER :: j,jj,k,l,m,info, setdim1, setdim, l2, k2
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
!
! Find indices of submatrix to be extracted and store them in set(1:setdim)
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel) .LT. 0) THEN
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
!
! Sort indices in set(1:setdim) in increasing order
!   (makes susequent processing faster)
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l)<set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
!
! Extract submatrix of (the symmetric part of) input matrix and store it in W
! Store in AGe the mean row and column sum of input matrix for related indices
! (this is also the row sum of A_G)
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=a1(idiag1(jj))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=a1(k)
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
    END DO
!
    DO j=1,SetDim
   !       Set sig(j) equal to minus the sum of connections that are
   !       external to the tentative aggregate; that is, add the internal
   !       offdiagonal entries to current sig(j)
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
   !       Heuristic if sig(j) is negative: set Sigma_G(j)=abs(sig(j))
   !         Hence need to correct the row sum of A_G stored in AGe
       IF (sig(j) < 0.0d0)  AGe(j)=AGe(j)+2*sig(j)
   !
   !       Set W = A_G ; offidagonal entries OK, need only to correct
   !       the diagonal
   !
       W(j,j)=W(j,j)-abs(sig(j))
   !
   !                    kap             1
   !       Now set W = ------ A_G  -  -----  M_G
   !                   kap-1          kap-1
   !
   !       Then W contains the matrix that need to be nonnegative definite to
   !       accept the aggregate according to the theory in [2],
   !       up to the rank one correction that is added below
   !
   !        Offdiagonal entrie of W alraedy OK [rule: offdiag(A_G)=ofdiag(MG)];
   !        for the diagonal ones, note that M_G(j,j)=A_G(j,j)+2*abs(sig(j)),
   !        hence the rule below with bndmum1m1=1/(kap-1)
   !
       tmp=2*abs(sig(j))
       W(j,j)=W(j,j)-bndmum1m1*tmp
   !
   !       Store M_G*1 = rowsum(A_G)+2*sig in v ; [ 1 = vector of all ones ];
   !       beta <- quadratic form (1^T * D_G * 1)
   !       alp is positive  if and only if AGe(j) is zero for all j;
   !              that is, if and only if A_G has 1 in its null space
       v(j)=tmp+AGe(j)
       IF (j .eq. 1) THEN
          beta=v(j)
          alp=abs(AGe(j))
       ELSE
          beta=beta+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
!
! Eventually, add to W the rank one correction
!                   1
!       -------------------------  M_G * 1 * 1^T * M_G
!       (kap-1)*(1^T * M_G * 1)
!
    beta=bndmum1m1/beta
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
!
! Now, the rule is to accept the tentative aggregate if and only if
! W(setdim,setdim) is nonnegative definite
! Two cases:
!     alp == 0 , then W is singular, hence we check the positive
!                definitenes of W(setdim-1,setdim-1)
!     alp > 0  , then W is normally not singular and we require thus that
!                it is positive definite
!
    IF (alp.LT.repsmach*beta) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
!
! Perform a Choleky factorization of W(setdim1,setdim1) and check
! for negative pivots; code expended for small dimensions (usual case)
!
    acc=.FALSE.
! acc has just be defined to FALSE; hence RETURN statement means rejection
!
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info .NE. 0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8) .LE. 0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7) .LE. 0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6) .LE. 0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5) .LE. 0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4) .LE. 0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3) .LE. 0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2) .LE. 0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1) .LE. 0.0d0) RETURN
10  CONTINUE
!
!   all test passed: accept the aggregate
    acc=.TRUE.
!
    RETURN
  END SUBROUTINE dagmg_checktentagg_SF
  END SUBROUTINE dagmg_findpairs_SF
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_GI(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
! Performs pairwise aggregation according to Algorithms 4.2 and 4.3 in [2,3],
! with some heuristic enhancements to cover non diagonally  dominant matrices.
!
! Version: initial pairwise aggregation (next levels) for general matrices [3].
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
!
!  INPUT: n,ja(*),ia(n+1),idiag(n): aggregated matrix from previous pass
!           (CSR format with partial ordering of each row: lower
!            triangular part, next diag., and eventually upper triagular part;
!            idiag point to the diagonal element)
!
!         si(n) : vector s or \tile{s} in Algorithm 4.2 or 4.3 of [2,3]
!
!         ind(n): on input, entries in ind should be nonnegative, and equal to
!             zero only at nodes that are to be kept outside the aggregation
!
!         nddl  : number of nodes i such that ind(i)=0
!
!         skipass: if true, pairwise aggregation is skipped; i.e.,
!               the coarsening only removes nodes i such that ind(i)==0
!
!         ipc: ipc(i) = 0 means that node i has to be transferred unaggregated
!                        to the coarse grid, and lcg(2,j) set to 0
!                        for the corresponding coarse grid node j
!
!         zerors: if true, diagonal entries are not taken into account;
!                 instead, the matrix is treated as it would have mean
!                 column and row sum equal to zero
!
!  OUTPUT: nc   : number of pairs formed
!
!          lcg(2,nc): list of aggregates;
!                     lcg(2,i)=0 : lcg(1,i) is a singleton that was
!                                  forced to be transferred unaggregated
!                     lcg(2,i)=-1: lcg(1,i) is a singleton because no valid
!                                  pair has been found
!                     lcg(2,i)>0 : (lcg(1,i),lcg(2,i)) form the ith pair
!
!          ind(n): ind(i)=j means that i belongs to aggregate j
!                  ind(i)=0 means that i does not belong to any pair
!                           and is listed to ldd
!                           (detected from ind(i)=0 on input)
!          ldd(nddl): list of nodes such that ind(i)=0 on input & output
!
!----------------
! Local variables
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       ijs=ijs+1
       !
       ! Node isel has been selected
       !
       ! First check if isel has to be kept outside aggregation
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       ! Check if isel has already been processed
       IF (ind(isel) .GE. 0) CYCLE
       !
       ! A new aggregate is formed that contains isel
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       ! Check if isel has to be transferred unaggregated
       IF (ipc(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ! Skip pairwise aggregation is skipass==true
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !  Try to form a pair with isel: follow list of neighbors
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          !  CYCLE means: reject this neighbor, proceed with next one
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          !  check if j is available to form a pair
          IF(ipc(j).EQ.0 .OR. ind(j).GE.0) CYCLE
          !   search for the corresponding entry in transposed matrix
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !   CYCLE instructions below: pair rejected because A_G is not
          !   nonnegative definite
          !
          !  Heuristic if sigj is negative: set Sigma_G(j)=abs(sigj) (j=1,2)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          !   main threshold test
          IF (valp > kaptg) CYCLE
          !
          !    A_G is nonneagtive definite and passed the corresponding
          !    "quality" threshold test: (isel,j) is
          !    an acceptable pair; check if it is the best one
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
       ! no valid pair found: isel left alone in aggregate nc
          lcg(2,nc) = -1
       ELSE
       ! pair nc is formed with isel and ipair
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
!
!   Restart aggregation with kaptg twice as large in case of
!   slow coarsening; allowed only once per level
!   (new parameters apply at further levels)
!
    IF (ifirst .AND. 3*nc.GE.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       ifirst=.FALSE.
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE dagmg_findpairs_GI
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_SI(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
! Performs pairwise aggregation according to Algorithms 4.2 and 4.3 in [2,3],
! with some heuristic enhancements to cover non diagonally  dominant matrices.
!
! Version: initial pairwise aggregation (next levels) for symm. matrices [2].
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
!
!  INPUT: n,ja(*),ia(n+1),idiag(n): aggregated matrix from previous pass
!           (CSR format with partial ordering of each row: lower
!            triangular part, next diag., and eventually upper triagular part;
!            idiag point to the diagonal element)
!
!         si(n) : vector s or \tile{s} in Algorithm 4.2 or 4.3 of [2,3]
!
!         ind(n): on input, entries in ind should be nonnegative, and equal to
!             zero only at nodes that are to be kept outside the aggregation
!
!         nddl  : number of nodes i such that ind(i)=0
!
!         skipass: if true, pairwise aggregation is skipped; i.e.,
!               the coarsening only removes nodes i such that ind(i)==0
!
!         ipc: ipc(i) = 0 means that node i has to be transferred unaggregated
!                        to the coarse grid, and lcg(2,j) set to 0
!                        for the corresponding coarse grid node j
!
!         zerors: if true, diagonal entries are not taken into account;
!                 instead, the matrix is treated as it would have mean
!                 column and row sum equal to zero
!
!  OUTPUT: nc   : number of pairs formed
!
!          lcg(2,nc): list of aggregates;
!                     lcg(2,i)=0 : lcg(1,i) is a singleton that was
!                                  forced to be transferred unaggregated
!                     lcg(2,i)=-1: lcg(1,i) is a singleton because no valid
!                                  pair has been found
!                     lcg(2,i)>0 : (lcg(1,i),lcg(2,i)) form the ith pair
!
!          ind(n): ind(i)=j means that i belongs to aggregate j
!                  ind(i)=0 means that i does not belong to any pair
!                           and is listed to ldd
!                           (detected from ind(i)=0 on input)
!          ldd(nddl): list of nodes such that ind(i)=0 on input & output
!
!----------------
! Local variables
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       ijs=ijs+1
       !
       ! Node isel has been selected
       !
       ! First check if isel has to be kept outside aggregation
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       ! Check if isel has already been processed
       IF (ind(isel) .GE. 0) CYCLE
       !
       ! A new aggregate is formed that contains isel
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       ! Check if isel has to be transferred unaggregated
       IF (ipc(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ! Skip pairwise aggregation is skipass==true
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !  Try to form a pair with isel: follow list of neighbors
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          !  CYCLE means: reject this neighbor, proceed with next one
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          !  check if j is available to form a pair
          IF(ipc(j).EQ.0 .OR. ind(j).GE.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !   CYCLE instructions below: pair rejected because A_G is not
          !   nonnegative definite
          !
          !  Heuristic if sigj is negative: set Sigma_G(j)=abs(sigj) (j=1,2)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          !   main threshold test
          IF (valp > kaptg) CYCLE
          !
          !    A_G is nonneagtive definite and passed the corresponding
          !    "quality" threshold test: (isel,j) is
          !    an acceptable pair; check if it is the best one
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
       ! no valid pair found: isel left alone in aggregate nc
          lcg(2,nc) = -1
       ELSE
       ! pair nc is formed with isel and ipair
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
!
!   Restart aggregation with kaptg twice as large in case of
!   slow coarsening; allowed only once per level
!   (new parameters apply at further levels)
!
    IF (ifirst .AND. 3*nc.GE.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       ifirst=.FALSE.
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE dagmg_findpairs_SI
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_GI1(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
! Performs pairwise aggregation according to Algorithms 4.2 and 4.3 in [2,3],
! with some heuristic enhancements to cover non diagonally  dominant matrices.
!
! Version: initial pairwise aggregation (top level) for general matrices [3].
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
!
!  INPUT: n,ja(*),ia(n+1),idiag(n): aggregated matrix from previous pass
!           (CSR format with partial ordering of each row: lower
!            triangular part, next diag., and eventually upper triagular part;
!            idiag point to the diagonal element)
!
!         si(n) : vector s or \tile{s} in Algorithm 4.2 or 4.3 of [2,3]
!
!         ind(n): on input, entries in ind should be nonnegative, and equal to
!             zero only at nodes that are to be kept outside the aggregation
!
!         nddl  : number of nodes i such that ind(i)=0
!
!         skipass: if true, pairwise aggregation is skipped; i.e.,
!               the coarsening only removes nodes i such that ind(i)==0
!
!         iperm(n),riperm(n): CMK and reverse CMK permutation, respectively
!           iperm(i) = 0 means that node i has to be transferred unaggregated
!                        to the coarse grid, and lcg(2,j) set to 0
!                        for the corresponding coarse grid node j
!
!         zerors: if true, diagonal entries are not taken into account;
!                 instead, the matrix is treated as it would have mean
!                 column and row sum equal to zero
!
!  OUTPUT: nc   : number of pairs formed
!
!          lcg(2,nc): list of aggregates;
!                     lcg(2,i)=0 : lcg(1,i) is a singleton that was
!                                  forced to be transferred unaggregated
!                     lcg(2,i)=-1: lcg(1,i) is a singleton because no valid
!                                  pair has been found
!                     lcg(2,i)>0 : (lcg(1,i),lcg(2,i)) form the ith pair
!
!          ind(n): ind(i)=j means that i belongs to aggregate j
!                  ind(i)=0 means that i does not belong to any pair
!                           and is listed to ldd
!                           (detected from ind(i)=0 on input)
!          ldd(nddl): list of nodes such that ind(i)=0 on input & output
!
!----------------
! Local variables
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       !
       ! Node isel has been selected
       !
       ! First check if isel has to be kept outside aggregation
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       ! Check if isel has already been processed
       IF (ind(isel) .GE. 0) CYCLE
       !
       ! A new aggregate is formed that contains isel
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       ! Check if isel has to be transferred unaggregated
       IF (iperm(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ! Skip pairwise aggregation is skipass==true
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !  Try to form a pair with isel: follow list of neighbors
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          !  CYCLE means: reject this neighbor, proceed with next one
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          !  check if j is available to form a pair
          IF(iperm(j).EQ.0 .OR. ind(j).GE.0) CYCLE
          !   search for the corresponding entry in transposed matrix
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk .NE. 0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !   CYCLE instructions below: pair rejected because A_G is not
          !   nonnegative definite
          !
          !  Heuristic if sigj is negative: set Sigma_G(j)=abs(sigj) (j=1,2)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          !   main threshold test
          IF (valp > kaptg) CYCLE
          !
          !    A_G is nonneagtive definite and passed the corresponding
          !    "quality" threshold test: (isel,j) is
          !    an acceptable pair; check if it is the best one
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
       ! no valid pair found: isel left alone in aggregate nc
          lcg(2,nc) = -1
       ELSE
       ! pair nc is formed with isel and ipair
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
!
!   Restart aggregation with kaptg twice as large in case of
!   slow coarsening; allowed only once per level
!   (new parameters apply at further levels)
!
    IF (ifirst .AND. 3*nc.GE.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       ifirst=.FALSE.
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE dagmg_findpairs_GI1
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_findpairs_SI1(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
! Performs pairwise aggregation according to Algorithms 4.2 and 4.3 in [2,3],
! with some heuristic enhancements to cover non diagonally  dominant matrices.
!
! Version: initial pairwise aggregation (top level) for symm. matrices [2].
!
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    REAL(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
!
!  INPUT: n,ja(*),ia(n+1),idiag(n): aggregated matrix from previous pass
!           (CSR format with partial ordering of each row: lower
!            triangular part, next diag., and eventually upper triagular part;
!            idiag point to the diagonal element)
!
!         si(n) : vector s or \tile{s} in Algorithm 4.2 or 4.3 of [2,3]
!
!         ind(n): on input, entries in ind should be nonnegative, and equal to
!             zero only at nodes that are to be kept outside the aggregation
!
!         nddl  : number of nodes i such that ind(i)=0
!
!         skipass: if true, pairwise aggregation is skipped; i.e.,
!               the coarsening only removes nodes i such that ind(i)==0
!
!         iperm(n),riperm(n): CMK and reverse CMK permutation, respectively
!           iperm(i) = 0 means that node i has to be transferred unaggregated
!                        to the coarse grid, and lcg(2,j) set to 0
!                        for the corresponding coarse grid node j
!
!         zerors: if true, diagonal entries are not taken into account;
!                 instead, the matrix is treated as it would have mean
!                 column and row sum equal to zero
!
!  OUTPUT: nc   : number of pairs formed
!
!          lcg(2,nc): list of aggregates;
!                     lcg(2,i)=0 : lcg(1,i) is a singleton that was
!                                  forced to be transferred unaggregated
!                     lcg(2,i)=-1: lcg(1,i) is a singleton because no valid
!                                  pair has been found
!                     lcg(2,i)>0 : (lcg(1,i),lcg(2,i)) form the ith pair
!
!          ind(n): ind(i)=j means that i belongs to aggregate j
!                  ind(i)=0 means that i does not belong to any pair
!                           and is listed to ldd
!                           (detected from ind(i)=0 on input)
!          ldd(nddl): list of nodes such that ind(i)=0 on input & output
!
!----------------
! Local variables
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       !
       ! Node isel has been selected
       !
       ! First check if isel has to be kept outside aggregation
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       ! Check if isel has already been processed
       IF (ind(isel) .GE. 0) CYCLE
       !
       ! A new aggregate is formed that contains isel
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       ! Check if isel has to be transferred unaggregated
       IF (iperm(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ! Skip pairwise aggregation is skipass==true
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !  Try to form a pair with isel: follow list of neighbors
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          !  CYCLE means: reject this neighbor, proceed with next one
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          !  check if j is available to form a pair
          IF(iperm(j).EQ.0 .OR. ind(j).GE.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !   CYCLE instructions below: pair rejected because A_G is not
          !   nonnegative definite
          !
          !  Heuristic if sigj is negative: set Sigma_G(j)=abs(sigj) (j=1,2)
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          !   main threshold test
          IF (valp > kaptg) CYCLE
          !
          !    A_G is nonneagtive definite and passed the corresponding
          !    "quality" threshold test: (isel,j) is
          !    an acceptable pair; check if it is the best one
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
       ! no valid pair found: isel left alone in aggregate nc
          lcg(2,nc) = -1
       ELSE
       ! pair nc is formed with isel and ipair
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
!
!   Restart aggregation with kaptg twice as large in case of
!   slow coarsening; allowed only once per level
!   (new parameters apply at further levels)
!
    IF (ifirst .AND. 3*nc.GE.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       ifirst=.FALSE.
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE dagmg_findpairs_SI1
!------------------------------------------------------------------
  SUBROUTINE dagmg_lcgmix(nc,m,lcg1,lcg,lcgn)
    INTEGER :: nc,m,lcg1(m,*),lcg(2,*),lcgn(2*m,*),i,l,l1,l2
    IF (m.eq.2) THEN
       DO i=1,nc
          IF(lcg(2,i) .EQ. 0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(4,i)=-1
          ELSE IF(lcg(2,i) .LT. 0) THEN
             IF (lcg1(2,lcg(1,i)) .LT. 0) THEN
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=-1
                lcgn(4,i)=-1
             ELSE
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=lcg1(2,lcg(1,i))
                lcgn(4,i)=-2
             END IF
          ELSE
             IF (lcg1(2,lcg(1,i)) .LT. 0) THEN
                IF (lcg1(2,lcg(2,i)) .LT. 0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-2
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(3,i)=lcg1(2,lcg(2,i))
                   lcgn(4,i)=-3
                END IF
             ELSE
                IF (lcg1(2,lcg(2,i)) .LT. 0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-3
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=lcg1(2,lcg(2,i))
                END IF
             END IF
          END IF
       END DO
    ELSE
       DO i=1,nc
          IF(lcg(2,i) .EQ. 0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(2*m,i)=-1
          ELSE
             lcgn(2,i)=-1
             l1=m
             IF (lcg1(m,lcg(1,i)).LT.0) l1=-lcg1(m,lcg(1,i))
             lcgn(1:l1,i)=lcg1(1:l1,lcg(1,i))
             IF(lcg(2,i) .LT. 0) THEN
                l=l1
             ELSE
                l2=m
                IF (lcg1(m,lcg(2,i)).LT.0) l2=-lcg1(m,lcg(2,i))
                lcgn(l1+1:l1+l2,i)=lcg1(1:l2,lcg(2,i))
                l=l1+l2
             END IF
             IF(l .LT. 2*m) lcgn(2*m,i)=-l
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE dagmg_lcgmix
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_setind(nc,ndd,ldd,lcg,m,ind)
    INTEGER :: nc,m,lcg(m,*),nll,ldd(ndd),ind(*),i,k,l
    DO i=1,ndd
       ind(ldd(i))=0
    END DO
    DO i=1,nc
       l=m
       IF (lcg(m,i) .LT. 0) l=-lcg(m,i)
       DO k=1,l
          ind(lcg(k,i))=i
       END DO
    END DO
    RETURN
  END SUBROUTINE dagmg_setind
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_setcg(n,a,ja,ia,idiag,si,ind,lcg           &
       ,nc,a2,ja2,ia2,idiag2,si2,ysi,maxdg,iw,w,iext,iext2)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),nc,lcg(2,nc)
    INTEGER :: ja2(*),ia2(n+1),idiag2(n),maxdg
    INTEGER, TARGET :: iw(2*nc)
    INTEGER, OPTIONAL :: iext(*),iext2(*)
    REAL(kind(0.0d0)) :: a(*),a2(*),w(nc),vald
    REAL(kind(0.0d0)) :: si(n),si2(*),sii
    LOGICAL :: ysi
    INTEGER :: nz,nzu,i,jj,jc,jcol,ki,kb,kf,jpos
    INTEGER, POINTER, DIMENSION(:) :: iw1, iw2
    !
    iw1 => iw(1:nc)
    iw2 => iw(nc+1:2*nc)
    !
    nz = 0
    iw1(1:nc)=0
    maxdg=0
    ia2(1)=1
    DO i = 1,nc
       sii=0.0d0
       vald=0.0d0
       nzu=0
       DO ki= 1,2
          jj = lcg(ki,i)
          IF (ki.EQ.1 .OR. jj.GT.0) THEN
             IF (ysi) sii=sii+si(jj)
             kf = ia(jj+1) - 1
             DO kb = ia(jj),kf
                jc = ja(kb)
                jcol = ind(jc)
                IF (jcol .GT. 0) THEN
                   IF (jcol .LT. i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nz = nz+1
                         ja2(nz) = jcol
                         iw1(jcol) = nz
                         a2(nz) = a(kb)
                      ELSE
                         a2(jpos) = a2(jpos) + a(kb)
                      ENDIF
                   ELSE IF (jcol .GT. i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nzu = nzu+1
                         iw2(nzu) = jcol
                         iw1(jcol) = nzu
                         w(nzu) = a(kb)
                      ELSE
                         w(jpos) = w(jpos) + a(kb)
                      ENDIF
                   ELSE
                      vald=vald+a(kb)
                      IF (ysi .AND. jc.NE.jj) sii=sii+a(kb)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       nz=nz+1
       a2(nz)=vald
       idiag2(i)=nz
       ja2(nz)=i
       a2(nz+1:nz+nzu)=w(1:nzu)
       ja2(nz+1:nz+nzu)=iw2(1:nzu)
       nz=nz+nzu
       maxdg=max(maxdg,nz-ia2(i))
       DO kb = ia2(i), nz
          iw1(ja2(kb))=0
       ENDDO
       IF (ysi) si2(i)=sii
       ia2(i+1)=nz+1
    ENDDO
    RETURN
  END SUBROUTINE dagmg_setcg
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of ROUTINES FOR SETUP
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! ROUTINES FOR DIRECT SOLUTION
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_LAPACK(N,a,ja,ia,f,ijb,flop)
    USE dagmg_mem
    IMPLICIT NONE
    INTEGER :: N,ia(n+1),ja(*),ijb
    REAL(kind(0.0d0)) :: a(*),f(n)
    REAL(kind(0.0d0)) :: flop
    !
    REAL(kind(0.0d0)), ALLOCATABLE, SAVE :: ac(:,:)
    INTEGER , ALLOCATABLE, SAVE :: ipiv(:)
    INTEGER, SAVE :: iflop
    INTEGER :: i,kk,iext
    INTEGER  :: ierr
    INTEGER , parameter :: IONE=1
    !
    ierr=0
    IF (ijb == -2) THEN
       !
       DEALLOCATE (ac,ipiv)
       !
    ELSE IF (ijb == 1) THEN
       !
       ALLOCATE (ac(n,n),ipiv(n))
       memi=memi+n
       memr=memr+n*n
       memax=MAX(memax,memr+memi*rlenilen)
       ac=0.0d0
       DO i=1,n
          DO kk=ia(i),ia(i+1)-1
             ac(i,ja(kk))=a(kk)
          END DO
       END DO
       CALL DGETRF(N,N,ac,N,ipiv,ierr)
       IF (ierr /= 0) THEN
          WRITE(iout, *) ' FATAL ERROR in GETRF: ierror=',ierr
          STOP
       END IF
       iflop=(2*n*n-n)
       flop=(2*1.0d0)/(3*1.0d0)*(dble(n)**3)
       !
    ELSE IF (ijb == 2) THEN
       !
       CALL DGETRS('N',N,IONE,ac,N,ipiv,f,N,ierr)
       IF (ierr /= 0) THEN
          WRITE(iout, *) ' FATAL ERROR in GETRS: ierror=',ierr
          STOP
       END IF
       flop=flop+dble(iflop)
       !
    END IF
    !
    RETURN
  END SUBROUTINE dagmg_LAPACK
!-----------------------------------------------------------------------
  SUBROUTINE dagmg_MUMPSseq(n,a,ja,ia,f,ijob,flop)
    USE dagmg_mem
    IMPLICIT NONE
!!!!!!!! BEGIN Include DAGMG_MUMPS definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      TYPE DMUMPS_ROOT_STRUC
        SEQUENCE
        INTEGER MBLOCK, NBLOCK, NPROW, NPCOL
        INTEGER MYROW, MYCOL
        INTEGER ROOT_SIZE, TOT_ROOT_SIZE
        INTEGER :: CNTXT_BLACS, truc
        INTEGER, DIMENSION(:), POINTER :: RG2L_ROW
        INTEGER, DIMENSION(:), POINTER :: RG2L_COL
        INTEGER , DIMENSION(:), POINTER :: IPIV
        INTEGER, DIMENSION( 9 ) :: DESCRIPTOR, DESCB
        LOGICAL yes, gridinit_done
        INTEGER LPIV, brol
!       Used to access Schur easily from root structure
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR_POINTER
        INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD, machin
!
!      Data for nullspace/QR
!
        DOUBLE PRECISION, DIMENSION(:), POINTER :: QR_TAU
        DOUBLE PRECISION     QR_RCOND
!
!      Givens rotations
!
        INTEGER MAXG, GIND
        DOUBLE PRECISION, DIMENSION(:),POINTER::GROW, GCOS, GSIN
!
!      RRRLU data
!
        INTEGER ELG_MAX,NULL_MAX
        INTEGER ELIND,EUIND,NLUPDATE,NUUPDATE
        INTEGER,DIMENSION(:),POINTER::PERM_ROW,PERM_COL
        INTEGER,DIMENSION(:),POINTER::ELROW, EUROW, PTREL, PTREU
        DOUBLE PRECISION, DIMENSION(:), POINTER :: ELELG, EUELG, DL
!
      END TYPE DMUMPS_ROOT_STRUC
      TYPE DMUMPS_STRUC
        SEQUENCE
!
! This structure contains all parameters
! for the interface to the user, plus internal
! information
!
! *****************
! INPUT PARAMETERS
! *****************
!    -----------------
!    MPI Communicator
!    -----------------
        INTEGER COMM
!    ------------------
!    Problem definition
!    ------------------
!    Solver (SYM=0 unsymmetric,SYM=1 symmetric Positive Definite,
!        SYM=2 general symmetric)
!    Type of parallelism (PAR=1 host working, PAR=0 host not working)
        INTEGER SYM, PAR
        INTEGER JOB
!    --------------------
!    Order of Input matrix
!    --------------------
        INTEGER N
!
!    ----------------------------------------
!    Assembled input matrix : User interface
!    ----------------------------------------
        INTEGER NZ
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A
        INTEGER, DIMENSION(:), POINTER :: IRN, JCN
        DOUBLE PRECISION, DIMENSION(:), POINTER :: COLSCA, ROWSCA, pad0
!
!       ------------------------------------
!       Case of distributed assembled matrix
!       matrix on entry:
!       ------------------------------------
        INTEGER NZ_loc, pad1
        INTEGER, DIMENSION(:), POINTER :: IRN_loc, JCN_loc
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A_loc, pad2
!
!    ----------------------------------------
!    Unassembled input matrix: User interface
!    ----------------------------------------
        INTEGER NELT, pad3
        INTEGER, DIMENSION(:), POINTER :: ELTPTR
        INTEGER, DIMENSION(:), POINTER :: ELTVAR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A_ELT, pad4
!
!    ---------------------------------------------
!    Symmetric permutation :
!               PERM_IN if given by user (optional)
!    ---------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: PERM_IN
!
!
! ******************
! INPUT/OUTPUT data
! ******************
!    --------------------------------------------------------
!    RHS / SOL_LOC
!    -------------
!       right-hand side and solution
!    -------------------------------------------------------
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS, REDRHS
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS_SPARSE
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SOL_LOC
        INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE
        INTEGER, DIMENSION(:), POINTER :: IRHS_PTR
        INTEGER, DIMENSION(:), POINTER :: ISOL_LOC
        INTEGER LRHS, NRHS, NZ_RHS, LSOL_LOC, LREDRHS
        INTEGER pad5
!    ----------------------------
!    Control parameters,
!    statistics and output data
!    ---------------------------
        INTEGER ICNTL(40)
        INTEGER INFO(40)
        INTEGER INFOG(40)
        DOUBLE PRECISION COST_SUBTREES
        DOUBLE PRECISION CNTL(15)
        DOUBLE PRECISION RINFO(20)
        DOUBLE PRECISION RINFOG(20)
!    ---------------------------------------------------------
!    Permutations computed during analysis:
!       SYM_PERM: Symmetric permutation
!       UNS_PERM: Column permutations (optionnal)
!    ---------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: SYM_PERM, UNS_PERM
!
!    -------------------------------------
!    Case of distributed matrix on entry:
!    DMUMPS potentially provides mapping
!    -------------------------------------
        INTEGER, DIMENSION(:), POINTER :: MAPPING
!
!    -------------------------------
!    Deficiency and null space basis
!    -------------------------------
        DOUBLE PRECISION, DIMENSION(:,:), POINTER :: NULL_SPACE
        INTEGER Deficiency, pad6
!    -----
!    Schur
!    -----
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK
        INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
        INTEGER SIZE_SCHUR
        INTEGER, DIMENSION(:), POINTER :: LISTVAR_SCHUR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR_CINTERFACE
!    --------------
!    Version number
!    --------------
        CHARACTER(LEN=16) VERSION_NUMBER
!    -----------
!    Out-of-core
!    -----------
        CHARACTER(LEN=256) :: OOC_TMPDIR
        CHARACTER(LEN=64) :: OOC_PREFIX
!    ------------------------------------------
!    To save the matrix in matrix market format
!    ------------------------------------------
        CHARACTER(LEN=256) WRITE_PROBLEM
!
!
! **********************
! INTERNAL Working data
! *********************
        INTEGER INST_Number
!       For MPI
        INTEGER COMM_NODES, MYID_NODES, COMM_LOAD
        INTEGER  MYID, NPROCS, NSLAVES
        INTEGER ASS_IRECV
        INTEGER, DIMENSION(:), POINTER :: POIDS
        INTEGER LBUFR
        INTEGER LBUFR_BYTES
        INTEGER, DIMENSION(:), POINTER ::  BUFR
!       For analysis/facto/solve phases
        INTEGER MAXIS1, pad7
        INTEGER KEEP(500)
        INTEGER*8 KEEP8(150)
!       IS is used for the factors + workspace for contrib. blocks
        INTEGER, DIMENSION(:), POINTER :: IS
!       is1 (maxis1) contains working arrays computed
!       and used only during analysis
        INTEGER, DIMENSION(:), POINTER :: IS1
!       The following data/arrays are computed during the analysis
!       phase and used during the factorization and solve phases.
        INTEGER LNA
        INTEGER NBSA
        INTEGER,POINTER,DIMENSION(:)::STEP, NE_STEPS, ND_STEPS
!  Info for pruning tree
        INTEGER,POINTER,DIMENSION(:)::Step2node
!  ---------------------
        INTEGER,POINTER,DIMENSION(:)::FRERE_STEPS, DAD_STEPS
        INTEGER,POINTER,DIMENSION(:)::FILS, PTRAR, FRTPTR, FRTELT
        INTEGER,POINTER,DIMENSION(:)::NA, PROCNODE_STEPS
!       The two pointer arrays computed in facto and used by the solve
!          (except the factors) are PTLUST_S and PTRFAC.
        INTEGER, DIMENSION(:), POINTER :: PTLUST_S
        INTEGER(8), DIMENSION(:), POINTER :: PTRFAC
!       main real working arrays for factorization/solve phases
        DOUBLE PRECISION, DIMENSION(:), POINTER :: S
!       Information on mapping
        INTEGER, DIMENSION(:), POINTER :: PROCNODE
!       Input matrix ready for numerical assembly
!           -arrowhead format in case of assembled matrix
!           -element format otherwise
        INTEGER, DIMENSION(:), POINTER :: INTARR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: DBLARR
!       Element entry: internal data
        INTEGER NELT_LOC, LELTVAR, NA_ELT, pad8
        INTEGER, DIMENSION(:), POINTER :: ELTPROC
!       Candidates and node partitionning
        INTEGER, DIMENSION(:,:), POINTER :: CANDIDATES
        INTEGER, DIMENSION(:),   POINTER :: ISTEP_TO_INIV2
        INTEGER, DIMENSION(:),   POINTER :: FUTURE_NIV2
        INTEGER, DIMENSION(:,:), POINTER :: TAB_POS_IN_PERE
        LOGICAL, DIMENSION(:), POINTER :: I_AM_CAND
!       For heterogeneous architecture
        INTEGER, DIMENSION(:), POINTER :: MEM_DIST
!       Compressed RHS
        INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHSCOMP
!       For C interface
!   Info on the subtrees to be used during factorization
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: MEM_SUBTREE
        INTEGER, DIMENSION(:),   POINTER :: MY_ROOT_SBTR
        INTEGER, DIMENSION(:),   POINTER :: MY_FIRST_LEAF
        INTEGER, DIMENSION(:),   POINTER :: MY_NB_LEAF
        INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: COST_TRAV
        INTEGER NBSA_LOCAL, zwave
        INTEGER(8) :: MAX_SURF_MASTER
        INTEGER :: LWK_USER, zozo
        DOUBLE PRECISION, DIMENSION(:), POINTER :: WK_USER
!    For simulating parallel out-of-core stack.
        DOUBLE PRECISION, DIMENSION(:),POINTER ::CB_SON_SIZE
!   Instance number used/managed by the C/F77 interface
        INTEGER INSTANCE_NUMBER
!    OOC management data that must persist from factorization to solve.
        INTEGER OOC_MAX_NB_NODES_FOR_ZONE
        INTEGER, DIMENSION(:,:),   POINTER :: OOC_INODE_SEQUENCE
        INTEGER(8),DIMENSION(:,:), POINTER :: OOC_SIZE_OF_BLOCK
        INTEGER*8, DIMENSION(:,:),   POINTER :: OOC_VADDR
        INTEGER,DIMENSION(:), POINTER :: OOC_TOTAL_NB_NODES
        INTEGER,DIMENSION(:), POINTER :: OOC_NB_FILES
        CHARACTER,DIMENSION(:,:), POINTER :: OOC_FILE_NAMES
        INTEGER,DIMENSION(:), POINTER :: OOC_FILE_NAME_LENGTH
!    Indices of nul pivots
        INTEGER,DIMENSION(:), POINTER :: PIVNUL_LIST
!    Internal control array
        DOUBLE PRECISION DKEEP(30)
!    Array needed to manage additionnal candidate processor
        INTEGER, DIMENSION(:,:), POINTER :: SUP_PROC
!   ------------------------
!   Root structure(internal)
!   ------------------------
        TYPE (DMUMPS_ROOT_STRUC) :: root
      END TYPE DMUMPS_STRUC
!!!!!!!! --END Include DAGMG_MUMPS definitions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER :: n,ia(n+1),ijob
    INTEGER, TARGET :: ja(*)
    REAL(kind(0.0d0)), TARGET :: a(*),f(n)
    REAL(kind(0.0d0)) :: flop
    !
    INTEGER, SAVE :: iflop
    TYPE(DMUMPS_STRUC), SAVE :: mumps_par
    INTEGER :: ierr, i, j, k
    !
    IF (ijob == -2) THEN
       !
       !  Destroy the instance (deallocate internal data structures)
       mumps_par%JOB = -2
       CALL DAGMG_MUMPS(mumps_par)
       !
    ELSE IF (ijob == 1) THEN
       !
       !  Define a communicator for the package.
       mumps_par%COMM = 0
       !  Initialize an instance of the package
       !  for L U factorization (sym = 0, with working host)
       mumps_par%JOB = -1
       mumps_par%SYM = 0
       mumps_par%PAR = 1
       CALL DAGMG_MUMPS(mumps_par)
       mumps_par%ICNTL(2)=-1
       mumps_par%ICNTL(3)=-1
       mumps_par%ICNTL(4)=0
       mumps_par%ICNTL(14)=80
       mumps_par%ICNTL(18)=3
       mumps_par%ICNTL(24)=1
       mumps_par%NZ_loc=ia(n+1)-1
       IF (transint) THEN
          ALLOCATE( mumps_par%JCN_loc(mumps_par%NZ_loc) )
          do i=1,n
             do j=ia(i),ia(i+1)-1
                mumps_par%JCN_loc(j)=i
             end do
          end do
          mumps_par%IRN_loc => ja(1:mumps_par%NZ_loc)
       ELSE
          ALLOCATE( mumps_par%IRN_loc(mumps_par%NZ_loc) )
          do i=1,n
             do j=ia(i),ia(i+1)-1
                mumps_par%IRN_loc(j)=i
             end do
          end do
          mumps_par%JCN_loc => ja(1:mumps_par%NZ_loc)
       END IF
       memi=memi+mumps_par%NZ_loc
       memax=MAX(memax,memr+memi*rlenilen)
       mumps_par%A_loc => a(1:mumps_par%NZ_loc)
       !
       mumps_par%N=n
       !  Call package for analysis+factorization
       mumps_par%JOB = 4
       CALL DAGMG_MUMPS(mumps_par)
       flop=mumps_par%RINFO(3)
       !iflop=2*mumps_par%INFO(9)-n
       IF (transint) THEN
          DEALLOCATE(mumps_par%JCN_loc)
          NULLIFY(mumps_par%IRN_loc,mumps_par%A_loc)
       ELSE
          DEALLOCATE(mumps_par%IRN_loc)
          NULLIFY(mumps_par%JCN_loc,mumps_par%A_loc)
       END IF
       memi=memi-mumps_par%NZ_loc
       !
    ELSE IF (ijob == 2) THEN
       !
       !  Call package for solution
       mumps_par%JOB = 3
       mumps_par%RHS => f(1:n)
       CALL DAGMG_MUMPS(mumps_par)
    END IF
    !
    RETURN
  END SUBROUTINE dagmg_MUMPSseq
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of ROUTINES FOR DIRECT SOLUTION
END MODULE dagmg_ALLROUTINES
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! MAIN DRIVER
!-----------------------------------------------------------------------
  SUBROUTINE dagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
    USE dagmg_mem
    USE dagmg_ALLROUTINES
    IMPLICIT NONE
    INTEGER    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
    REAL (kind(0.0d0)) :: a(*),f(n),x(n)
    REAL (kind(0.0d0)) :: tol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Arguments
!  =========
!
!  N       (input) INTEGER.
!          The dimension of the matrix.
!
!  A       (input/output) REAL (kind(0.0d0)). Numerical values of the matrix.
!  IA      (input/output) INTEGER. Pointers for every row.
!  JA      (input/output) INTEGER. Column indices.
!
!              AGMG ASSUMES THAT ALL DIAGONAL ENTRIES ARE POSITIVE
!
!          Detailed description of the matrix format
!
!              On input, IA(I), I=1,...,N, refers to the physical start
!              of row I. That is, the entries of row I are located
!              in A(K), where K=IA(I),...,IA(I+1)-1. JA(K) carries the
!              associated column indices. IA(N+1) must be defined in such
!              a way that the above rule also works for I=N (that is,
!              the last valid entry in arrays A,JA should correspond to
!              index K=IA(N+1)-1). According what is written
!              above, AGMG assumes that some of these JA(K) (for
!              IA(I)<= K < IA(I+1) ) is equal to I with corresponding
!              A(K) carrying the value of the diagonal element,
!              which must be positive.
!
!              A,IA,JA are "output" parameters because on exit the
!              entries of each row may occur in a different order (The
!              matrix is mathematically the same, but stored in
!              different way).
!
!  F       (input/output) REAL (kind(0.0d0)).
!          On input, the right hand side vector f.
!          Overwritten on output.
!          Significant only if IJOB==0, 2, 3, 10, 12, 100, 102, 110, 112
!
!  X       (input/output) REAL (kind(0.0d0)).
!          On input and if IJOB== 10, 12, 110, 112: initial guess
!             (for other values of IJOB, the default is used: the zero vector).
!          On output, the computed solution.
!
! IJOB     (input) INTEGER. Tells AGMG what has to be done.
!          0: performs setup + solve + memory release, no initial guess
!         10: performs setup + solve + memory release, initial guess in x(1:n)
!          1: performs setup only
!             (preprocessing: prepares all parameters for subsequent solves)
!          2: solves only (based on previous setup), no initial guess
!         12: solves only (based on previous setup), initial guess in x(1:n)
!          3: the vector returned in x(1:n) is not the solution of the linear
!                 system, but the result of the action of the multigrid
!                 preconditioner on the right hand side in f(1:n)
!         -1: erases the setup and releases internal memory
!
!   IJOB == 100,110,101,102,112: same as, respectively, IJOB==0,10,1,2,12
!       but, use the TRANSPOSE of the input matrix in A, JA, IA.
!
!   !!! IJOB==2,3,12,102,112 require that one has previously called AGMG
!       with IJOB==1 or IJOB==101
!
!   !!! (change with respect to versions 2.x) !!!
!       The preconditioner defined when calling AGMG
!         with IJOB==1 or IJOB==101 is entirely kept in internal memory.
!       Hence the arrays A, JA and IA are not accessed upon subsequent calls
!         with IJOB==3.
!       Upon subsequent calls with IJOB==2,12,102,112, a matrix needs to
!            be supplied in arrays A, JA, IA, but it will be used to
!            perform matrix vector product within the main iterative
!            solution process (and only for this).
!            Hence the system is solved with this matrix which
!            may differ from the matrix in A, JA, IA that was supplied
!            upon the previous call with IJOB==1 or IJOB==101;
!            then AGMG attempts to solve a linear system with the "new"
!            matrix (supplied when IJOB==2,12,102 or 112) using the
!            preconditioner set up for the "old" one (supplied when
!            IJOB==1 or 101).
!         The same remarks apply to IJOB >= 100 or not: the value IJOB==1
!            or 101 determines whether the preconditioner set up and stored
!            in internal memory is based on the matrix or its transpose;
!            the value IJOB==2,12 or 102,112 is used to determine whether
!            the linear system to be solved is with the matrix or its
!            transpose, independently of the set up.
!            Hence one may set up a preconditioner for a matrix and use it
!            for its transpose.
!       These functionalities (set up a preconditioner and use it for another
!            matrix) are provided for the sake of generality but should be
!            used with care; in general, set up is fast with AGMG and hence
!            it is recommended to rerun it even if the matrix changes only
!            slightly.
!
! IPRINT   (input) INTEGER.
!              Indicates the unit number where information is to be printed
!              (N.B.: 5 is converted to 6). If nonpositive, only error
!              messages are printed on standard output.
!
! NREST    (input) INTEGER.
!             Restart parameter for GCR (an implementation of GMRES)
!             Nonpositive values are converted to NREST=10 (default)
!
! !!  If NREST==1, Flexible CG is used instead of GCR (when IJOB==0,10,2,
!             12,100,110,102,112) and also (IJOB==0,1,100,101) performs some
!             simplification during the set up based on the assumption
!             that the matrix supplied in A, JA, IA is symmetric (there is
!             then no more difference between IJOB==1 and IJOB==101).
!
! !!!  NREST=1 Should be used if and only if the matrix is really SYMMETRIC
! !!!         (and positive definite).
!
!  ITER    (input/output) INTEGER.
!          On input, the maximum number of iterations. Should be positive.
!          On output, actual number of iterations.
!            If this number of iteration was insufficient to meet convergence
!            criterion, ITER will be returned negative and equal to the
!            opposite of the number of iterations performed.
!          Significant only if IJOB==0, 2, 10, 12, 100, 102, 110, 112
!
!  TOL     (input) REAL (kind(0.0d0)).
!          The tolerance on residual norm. Iterations are stopped whenever
!               || A*x-f || <= TOL* || f ||
!          Should be positive and less than 1.0
!          Significant only if IJOB==0, 2, 10, 12, 100, 102, 110, 112
!
!!!!! Remark !!!! Except insufficient number of iterations to achieve
!                 convergence (characterized by a negative value returned
!                 in ITER), all other detected errors are fatal and lead
!                 to a STOP statement.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Local variables
!
    INTEGER, SAVE :: nza
    REAL(kind(0.0d0)), SAVE :: cputm=0.0d0,eltm=0.0d0,flop=0.0d0
    LOGICAL, SAVE :: preprocessed=.FALSE.,solve=.FALSE.
    INTEGER :: i,init,ijb,k
    REAL(kind(0.0d0)) :: cputmp,eltmp,memh
    REAL(kind(0.0d0)) :: resid
    INTEGER, POINTER, DIMENSION(:) :: iw
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: ascal
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: w
    REAL(kind(0.0d0)) :: adum(1),fdum(1)
    INTEGER :: listrank(1),ifirstlistrank
!
    wfo=.TRUE.
    iout=iprint
    IF (iprint <= 0) THEN
       iout=6
       wfo=.FALSE.
    ELSE IF (iprint == 5) THEN
       iout=6
    END IF
    ijb=mod(ijob,100)
!  transpose on input if ijob >= 100
    trans=ijob.GE.100
!  matrix is assume symmetric positive definite if and only if nrst==1
    spd=nrest.EQ.1
!  transposition disabled at inner levels if matrix supposed spd
    transint=trans.AND.(.NOT.spd)
!
    wff=wfo.AND.(IRANK<=0)
!
    IF (MOD(ijb,10) >= 2 .AND. .NOT.preprocessed) THEN
       WRITE (iout,1001) IRANK,ijob
       STOP
    END IF
!
    IF (ijb < 0 .AND. solve) GOTO 450
    IF (ijb < 0) GOTO 500
    IF (MOD(ijb,10) >= 2) GOTO 300
    IF (preprocessed) THEN
       CALL dagmg_relmem
       IF (.NOT.allzero)                        &
               CALL dagmg_MUMPSseq(nn(nlev),adum,ja,ia,fdum,-2,flop)
       preprocessed=.FALSE.
       solve=.FALSE.
       eltm=0.0d0
       cputm=0.0d0
       flop=0.0d0
       kstat=0
    END IF
    CALL dagmg_mestime(-1,0.0d0,0.0d0)
!
!       Initial settings
!
    nza=ia(n+1)-ia(1)
    IF (HUGE(n) > 1.0e10) THEN
       rlenilen=dble(8)/dble(real_len)
    ELSE
       rlenilen=dble(4)/dble(real_len)
    END IF
!
!
!-----  PREPROCESSING
!
    nlev=0
    imult=1
    gcrin=.FALSE.
    IF (wfo) THEN
       WRITE(iout,900) IRANK
    END IF
!
!   Partial ordering of the rows
!
    ALLOCATE(dt(1)%idiag(n+1),w(n),iw(n))
    CALL dagmg_partroword(n,a,ja,ia,dt(1)%idiag,w,iw)
    DEALLOCATE(w,iw)
    !
    ! Setup: define hierarchy of grids and matrices
    !
    CALL dagmg_setupL1(n,a,ja,ia,listrank,ifirstlistrank)
    !
    ! Matrix is defined at each level: now perform smoother set up
    !
    CALL dagmg_smoothsetup
    preprocessed=.TRUE.
    memh=memr+memi*rlenilen
!
    CALL dagmg_mestime(1,cputmp,eltmp)
    IF(wfo)THEN
       IF (nlev > 1) THEN
          WRITE(iout,960) IRANK, memax/nza, real_len    &
                        , memax*real_len/(2**20)
          IF(MOD(ijb,10) == 1) THEN
             WRITE(iout,961) IRANK, memh/nza, real_len  &
                           , memh*real_len/(2**20)
          END IF
       END IF
   !CPU_TIME: next line may be uncommented if implemented
   !          (see subroutine mestime)
   !              WRITE(iout,996) IRANK,cputmp
       WRITE(iout,997) IRANK,eltmp
       WRITE(iout,'()')
    END IF
!
    IF (MOD(ijb,10) == 1) RETURN
!
!-----  SOLUTION PROCESS
!
300 CONTINUE
    CALL dagmg_mestime(-2,0.0d0,0.0d0)
    resid=tol
    nrst=nrest
    init=MAX(IJB/10,0)
    IF (nrst <= 0) nrst=10 ! Default restart paramater for GMRESR
   !
   !-     SINGLE APPLICATION OF MG PRECONDITIONER
    IF (MOD(ijb,10) >= 3) THEN
       IF (wfo) THEN
          WRITE(iout,901) IRANK
       END IF
       CALL dagmg_applyprec(n,f,x,a,ja,ia,flop)
       CALL dagmg_mestime(2,cputmp,eltmp)
       cputm=cputm+cputmp
       eltm=eltm+eltmp
       solve=.TRUE.
       RETURN
   !
   !-     END OF MG PREC
   !
   !-     GCR SOLUTION
    ELSE IF (nrst > 1) THEN
   !
       CALL dagmg_GCR(n,f,x,iter,resid,a,ja,ia,init,flop)
   !
   !-     END OF GCR SOLVE
   !
   !-     FLEXIBLE CG SOLUTION
    ELSE
   !
       CALL dagmg_FlexCG(n,f,x,iter,resid,a,ja,ia,init,flop)
   !
   !-     END OF FCG SOLVE
   !
    END IF
!
    CALL dagmg_mestime(2,cputm,eltm)
    solve=.FALSE.
    GOTO 455
450 CONTINUE
    IF (wfo) THEN
       WRITE(iout,'()')
       WRITE(iout,990) IRANK
       WRITE(iout,'()')
    END IF
455 CONTINUE
    IF (wfo) THEN
       IF (wff .AND. iter.NE.0) THEN
          DO i=2,nlev-1
             WRITE(iout,955) i,kstat(2,i-1),kstat(2,i),         &
                  dble(kstat(2,i))/dble(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       WRITE(iout,954) IRANK,flop/dble(11*n+2*nza)
       WRITE(iout,962) IRANK,(memh+mritr)/nza, real_len   &
                     , (memh+mritr)*real_len/(2**20)
   !CPU_TIME: next line may be uncommented if implemented
   !          (see subroutine mestime)
   !            WRITE(iout,998) IRANK,cputm
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
    END IF
       solve=.FALSE.
       eltm=0.0d0
       cputm=0.0d0
       flop=0.0d0
       kstat=0
!
    IF (MOD(ijb,10) > 0) RETURN
!
!----  RELEASE MEMORY
!
500 CONTINUE
    CALL dagmg_relmem
    IF (.NOT.allzero)                              &
            CALL dagmg_MUMPSseq(nn(nlev),adum,ja,ia,fdum,-2,flop)
    preprocessed=.FALSE.
    solve=.FALSE.
    eltm=0.0d0
    cputm=0.0d0
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE (iout,902) IRANK
    END IF
!
    CONTINUE
!
    RETURN
900 FORMAT(i3,'*ENTERING AGMG **********************************',&
         '***************************')
901 FORMAT(i3,'*ONE APPLICATION OF AGMG PRECONDITIONER')
902 FORMAT(i3,'*LEAVING AGMG * (MEMORY RELEASED) ***************',&
         '***************************')
903 FORMAT(i3,'*SYSTEM TOO SMALL: use direct solver')
954 FORMAT(i3,'*','  Equiv. number of CG iter.:',f9.2)
955 FORMAT('****     level',i2,'   #call=',i6,'   #cycle=',i6,    &
         '   mean=',f7.2,'    max=',i3)
960 FORMAT(  i3,'*','         memory used (peak):',f9.2,        &
         ' real(',i1,') words per nnz (',f8.2,' Mb)')
961 FORMAT(  i3,'*','     memory still allocated:',f9.2,        &
         ' real(',i1,') words per nnz (',f8.2,' Mb)')
962 FORMAT(  i3,'*','   memory used for solution:',f9.2,        &
         ' real(',i1,') words per nnz (',f8.2,' Mb)')
990 FORMAT(i3,'*GLOBAL STATISTICS for preconditioner application:')
996 FORMAT(i3,'*','           Setup time (CPU):   ',1PE10.2,     &
         ' seconds')
997 FORMAT(i3,'*','       Setup time (Elapsed):   ',1PE10.2,     &
         ' seconds')
998 FORMAT(i3,'*','        Solution time (CPU):   ',1PE10.2,     &
         ' seconds')
999 FORMAT(i3,'*','    Solution time (Elapsed):   ',1PE10.2,     &
         ' seconds')
1001 FORMAT(i3,'*',' FATAL ERROR: setup not done: ijob=',i3, &
         ' is not allowed')
1002 FORMAT(i3,'*',' FATAL ERROR: ijob=',i3, &
         ' (i.e. >= 100: work with transpose) is not allowed in the || case')
  END SUBROUTINE dagmg
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of MAIN DRIVER
!-----------------------------------------------------------------------
!------------------ END of source file ---------------------------------
