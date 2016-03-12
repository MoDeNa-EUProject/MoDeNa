!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! SOLVING A LINEAR MATRIX SYSTEM  AX = B  with Gauss-Jordan method using full
! pivoting at each step. During the process, original A and B matrices are
! destroyed to spare storage location.

! INPUTS:    A   MATRIX N*N
!            B   MATRIX N*M

! OUTPUS:    A   INVERSE OF A N*N
!            DET  DETERMINANT OF A
!            B   SOLUTION MATRIX N*M

! NOTA - If M=0 inversion of A matrix only.
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine MATINV ( N, M, AA, BB, DET )

  implicit none
  !-----------------------------------------------------------------------------
  integer, intent(in)                    :: N
  integer, intent(in)                    :: M
  real, dimension(N,N), intent(in out)   :: AA
  real, dimension(N,M), intent(in out)   :: BB
  real, intent(out)                      :: DET

  !-----------------------------------------------------------------------------
  integer                                :: i, j, k, ik, jk
  real, PARAMETER                        :: EPSMACH = 2.E-16
  real                                   :: PV, PAV, TT    
  integer, allocatable                   :: PC(:), PL(:)
  real, allocatable                      :: CS(:)
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! Initializations
  !-----------------------------------------------------------------------------
  allocate( PC(1:N) )
  allocate( PL(1:N) )
  allocate( CS(1:N) )
  DET = 1.0
  do I = 1, N
     PC(I) = 0
     PL(I) = 0
     CS(I) = 0.0
  end do

  !-----------------------------------------------------------------------------
  ! main loop
  !-----------------------------------------------------------------------------
  do K = 1, N
     !Searching greatest pivot
     PV = AA(K,K)
     IK = K
     JK = K
     PAV = ABS(PV)
     do I = K, N
        do J = K, N
           if ( ABS(AA(I,J)) > PAV ) then
              PV = AA(I,J)
              PAV = ABS(PV)
              IK = I
              JK = J
           end if
        end do
     end do

     !--------------------------------------------------------------------------
     ! Search terminated, the pivot is in location I=IK, J=JK.
     ! Memorizing pivot location:
     !--------------------------------------------------------------------------
     PC(K) = JK
     PL(K) = IK

     !--------------------------------------------------------------------------
     ! Determinant DET is actualised
     ! If DET=0, ERROR MESSAGE and STOP
     ! Machine dependant EPSMACH equals here 1.DE-20
     !--------------------------------------------------------------------------

     if ( IK /= K ) DET = -DET
     if ( JK /= K ) DET = -DET
     DET = DET * PV
     if ( ABS(DET) < EPSMACH ) then
        !Error message and Stop
        write (*,*) ' DETERMINANT EQUALS ZERO, NO SOLUTION!'
        STOP
     end if

     !--------------------------------------------------------------------------
     ! POSITIONNING PIVOT IN K,K
     !--------------------------------------------------------------------------
     if ( IK /= K ) then
        do I=1,N
           !EXCHANGE LINES IK and K
           TT=AA(IK,I)
           AA(IK,I)=AA(K,I)
           AA(K,I)=TT
        end do
     end if
     if (M /= 0) then
        do I = 1, M
           TT=BB(IK,I)
           BB(IK,I)=BB(K,I)
           BB(K,I)=TT
        end do
     end if

     !--------------------------------------------------------------------------
     ! Pivot is at correct line
     !--------------------------------------------------------------------------
     if (JK /= K) then
        do I = 1, N
           !Exchange columns JK and K of matrix AA
           TT = AA(I,JK)
           AA(I,JK) = AA(I,K)
           AA(I,K) = TT
        end do
     end if

     !--------------------------------------------------------------------------
     ! Pivot is at correct column and located in K,K
     !--------------------------------------------------------------------------
     ! Store column K in vector CS
     ! then set column K to zero
     !--------------------------------------------------------------------------
     do I = 1, N
        CS(I) = AA(I,K)
        AA(I,K) = 0.0
     end do

     CS(K) = 0.0
     AA(K,K) = 1.0
     !--------------------------------------------------------------------------
     ! Modify line K
     !--------------------------------------------------------------------------
     if (ABS(PV) < EPSMACH) then
        WRITE(*,*) '  PIVOT TOO SMALL - STOP'
        STOP
     end if
     do I = 1, N
        AA(K,I) = AA(K,I) / PV
     end do
     if (M /= 0) then
        do I = 1, M
           BB(K,I) = BB(K,I) / PV
        end do
     end if
     !--------------------------------------------------------------------------
     ! Modify other lines of matrix AA
     !--------------------------------------------------------------------------
     do J = 1, N
        if ( J == K ) CONTINUE
        do I = 1, N
           !Modify line J of matrix AA
           AA(J,I) = AA(J,I) - CS(J) * AA(K,I)
        end do
        if ( M /= 0 ) then
           do I = 1, M
              BB(J,I) = BB(J,I) - CS(J) * BB(K,I)
           end do
        end if
     end do
     !Line K is ready
  end do
  !End of K loop

  !-----------------------------------------------------------------------------
  ! The matrix AA is inverted - Rearrange AA
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !Exchange lines
  !-----------------------------------------------------------------------------
  do I = N, 1, -1
     IK = PC(I)
     if ( IK == I ) CONTINUE
     !EXCHANGE LINES I AND PC(I) OF AA
     do J = 1, N
        TT = AA(I,J)
        AA(I,J) = AA(IK,J)
        AA(IK,J) = TT
     end do
     if ( M /= 0 ) then
        do J = 1, M
           TT = BB(I,J)
           BB(I,J) = BB(IK,J)
           BB(IK,J) = TT
        end do
     end if
     !NO MORE EXCHANGE NEEDED
     !GO TO NEXT LINE
  end do


  !-----------------------------------------------------------------------------
  ! EXCHANGE COLUMNS
  !-----------------------------------------------------------------------------
  do J = N, 1, -1
     JK = PL(J)
     if ( JK == J ) CONTINUE
     !EXCHANGE COLUMNS J AND PL(J) OF AA
     do I = 1, N
        TT = AA(I,J)
        AA(I,J) = AA(I,JK)
        AA(I,JK) = TT
     end do
     !NO MORE EXCHANGE NEEDED
     !GO TO NEXT COLUMN
  end do

  deallocate( PC, PL, CS )

end subroutine MATINV

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine jacobi_eigenvalue
!
! This function computes the eigenvalues and eigenvectors of a real symmetric
! matrix, using Rutishauser's modfications of the classical Jacobi rotation
! method with threshold pivoting.
!
! input:
!    n        order of the matrix
!    a(n,n)   matrix, which must be square, real, and symmetric
!    it_max   maximum number of iterations
!
! output:
!    v(n,n)   matrix of eigenvectors
!    d(n)     eigenvalues, in descending order
!    it_num   total number of iterations
!    rot_num  total number of rotations
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )

  implicit none
  !-----------------------------------------------------------------------------
  integer, intent(in)                    :: n
  real, dimension(n,n), intent(in out)   :: a
  integer, intent(in)                    :: it_max
  real, dimension(n,n), intent(out)      :: v
  real, dimension(n), intent(out)        :: d
  integer, intent(out)                   :: it_num
  integer, intent(out)                   :: rot_num
  !-----------------------------------------------------------------------------

  integer                                :: i, j, k, l, m, p, q
  real                                   :: c
  real                                   :: g
  real                                   :: gapq
  real                                   :: h
  real                                   :: s
  real                                   :: t
  real                                   :: tau
  real                                   :: term
  real                                   :: termp
  real                                   :: termq
  real                                   :: theta
  real                                   :: thresh
  real, dimension(n)                     :: bw
  real, dimension(n)                     :: w, zw
  !-----------------------------------------------------------------------------

  do j = 1, n
     do i = 1, n
        v(i,j) = 0.0
     end do
     v(j,j) = 1.0
  end do

  do i = 1, n
     d(i) = a(i,i)
  end do

  bw(1:n) = d(1:n)
  zw(1:n) = 0.0
  it_num = 0
  rot_num = 0

  do while ( it_num < it_max )

     it_num = it_num + 1
     !--------------------------------------------------------------------------
     !  The convergence threshold is based on the size of the elements in
     !  the strict upper triangle of the matrix.
     !--------------------------------------------------------------------------
     thresh = 0.0
     do j = 1, n
        do i = 1, j - 1
           thresh = thresh + a(i,j) ** 2
        end do
     end do

     thresh = sqrt ( thresh ) / real ( 4 * n, kind = 8 )

     if ( thresh == 0.0 ) then
        exit 
     end if

     do p = 1, n
        do q = p + 1, n

           gapq = 10.0 * abs ( a(p,q) )
           termp = gapq + abs ( d(p) )
           termq = gapq + abs ( d(q) )
           !--------------------------------------------------------------------
           !  Annihilate tiny offdiagonal elements.
           !--------------------------------------------------------------------
           if ( 4 < it_num .and. termp == abs ( d(p) ) .and. termq == abs ( d(q) ) ) then

              a(p,q) = 0.0
              !-----------------------------------------------------------------
              !  Otherwise, apply a rotation.
              !-----------------------------------------------------------------
           else if ( thresh <= abs ( a(p,q) ) ) then

              h = d(q) - d(p)
              term = abs ( h ) + gapq

              if ( term == abs ( h ) ) then
                 t = a(p,q) / h
              else
                 theta = 0.5 * h / a(p,q)
                 t = 1.0 / ( abs ( theta ) + sqrt ( 1.0 + theta * theta ) )
                 if ( theta < 0.0 ) then 
                    t = - t
                 end if
              end if

              c = 1.0 / sqrt ( 1.0 + t * t )
              s = t * c
              tau = s / ( 1.0 + c )
              h = t * a(p,q)
              !-----------------------------------------------------------------
              !  Accumulate corrections to diagonal elements.
              !-----------------------------------------------------------------
              zw(p) = zw(p) - h                  
              zw(q) = zw(q) + h
              d(p) = d(p) - h
              d(q) = d(q) + h

              a(p,q) = 0.0
              !-----------------------------------------------------------------
              !  Rotate, using information from the upper triangle of A only.
              !-----------------------------------------------------------------
              do j = 1, p - 1
                 g = a(j,p)
                 h = a(j,q)
                 a(j,p) = g - s * ( h + g * tau )
                 a(j,q) = h + s * ( g - h * tau )
              end do

              do j = p + 1, q - 1
                 g = a(p,j)
                 h = a(j,q)
                 a(p,j) = g - s * ( h + g * tau )
                 a(j,q) = h + s * ( g - h * tau )
              end do

              do j = q + 1, n
                 g = a(p,j)
                 h = a(q,j)
                 a(p,j) = g - s * ( h + g * tau )
                 a(q,j) = h + s * ( g - h * tau )
              end do
              !-----------------------------------------------------------------
              !  Accumulate information in the eigenvector matrix.
              !-----------------------------------------------------------------
              do j = 1, n
                 g = v(j,p)
                 h = v(j,q)
                 v(j,p) = g - s * ( h + g * tau )
                 v(j,q) = h + s * ( g - h * tau )
              end do

              rot_num = rot_num + 1

           end if

        end do
     end do

     bw(1:n) = bw(1:n) + zw(1:n)
     d(1:n) = bw(1:n)
     zw(1:n) = 0.0

  end do
  !-----------------------------------------------------------------------------
  !  Restore upper triangle of input matrix.
  !-----------------------------------------------------------------------------
  do j = 1, n
     do i = 1, j - 1
        a(i,j) = a(j,i)
     end do
  end do
  !-----------------------------------------------------------------------------
  !  Ascending sort the eigenvalues and eigenvectors.
  !-----------------------------------------------------------------------------
  do k = 1, n - 1

     m = k

     do l = k + 1, n
        if ( d(l) < d(m) ) then
           m = l
        end if
     end do

     if ( m /= k ) then

        t    = d(m)
        d(m) = d(k)
        d(k) = t

        w(1:n)   = v(1:n,m)
        v(1:n,m) = v(1:n,k)
        v(1:n,k) = w(1:n)

     end if

  end do

end subroutine jacobi_eigenvalue


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! subroutine eigenvalues ( n, a, x )
!
! evaluates eigenvalues and eigenvectors of a real symmetric matrix
!    a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices, Alex G. (December 2009)
!
! input:
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output:
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

!!$subroutine eigenvalues ( n, a, x )
!!$
!!$  implicit none
!!$  !-----------------------------------------------------------------------------
!!$  integer, intent(in)                    :: n
!!$  real, dimension(n,n), intent(in out)   :: a
!!$  real, dimension(n,n), intent(out)      :: x
!!$  !-----------------------------------------------------------------------------
!!$  integer                                :: i, j, k
!!$  real                                   :: b2, bar
!!$  real                                   :: beta, coeff, c, s, cs, sc
!!$  real, parameter                        :: abserr = 1.E-08
!!$  !-----------------------------------------------------------------------------
!!$  write (*,*) 'enter eigenvalues'
!!$  !-----------------------------------------------------------------------------
!!$  ! initialize x(i,j)=0, x(i,i)=1
!!$  !-----------------------------------------------------------------------------
!!$  x = 0.0
!!$  do i = 1, n
!!$     x(i,i) = 1.0
!!$  end do
!!$
!!$  !-----------------------------------------------------------------------------
!!$  ! find the sum of all off-diagonal elements (squared)
!!$  !-----------------------------------------------------------------------------
!!$  b2 = 0.0
!!$  do i = 1, n
!!$     do j = 1, n
!!$        if ( i /= j ) b2 = b2 + a(i,j)**2
!!$     end do
!!$  end do
!!$
!!$  if ( b2 <= abserr ) return
!!$
!!$  !-----------------------------------------------------------------------------
!!$  ! average for off-diagonal elements /2
!!$  !-----------------------------------------------------------------------------
!!$  bar = 0.5 * b2 / float(n*n)
!!$
!!$  do while ( b2 > abserr )
!!$     do i = 1, n-1
!!$        do j = i+1, n
!!$           if ( a(j,i)**2 <= bar ) cycle  ! do not touch small elements
!!$           b2 = b2 - 2.0*a(j,i)**2
!!$           bar = 0.5 * b2 / float(n*n)
!!$           !--------------------------------------------------------------------
!!$           ! calculate coefficient c and s for Givens matrix
!!$           !--------------------------------------------------------------------
!!$           beta = ( a(j,j) - a(i,i) ) / ( 2.0 * a(j,i) )
!!$           coeff = 0.5 * beta / sqrt( 1.0+beta**2 )
!!$           s = sqrt(max(0.5+coeff,0.0))
!!$           c = sqrt(max(0.5-coeff,0.0))
!!$           !--------------------------------------------------------------------
!!$           ! recalculate rows i and j
!!$           !--------------------------------------------------------------------
!!$           do k = 1, n
!!$              cs =  c * a(i,k) + s * a(j,k)
!!$              sc = -s * a(i,k) + c * a(j,k)
!!$              a(i,k) = cs
!!$              a(j,k) = sc
!!$           end do
!!$           !--------------------------------------------------------------------
!!$           ! new matrix a_{k+1} from a_{k}, and eigenvectors 
!!$           !--------------------------------------------------------------------
!!$           do k = 1, n
!!$              cs =  c * a(k,i) + s * a(k,j)
!!$              sc = -s * a(k,i) + c * a(k,j)
!!$              a(k,i) = cs
!!$              a(k,j) = sc
!!$              cs =  c * x(k,i) + s * x(k,j)
!!$              sc = -s * x(k,i) + c * x(k,j)
!!$              x(k,i) = cs
!!$              x(k,j) = sc
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  write (*,*) 'exit eigenvalues'
!!$end subroutine eigenvalues



