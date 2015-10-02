subroutine dgtsl ( n, c, d, e, b, info )

!*****************************************************************************80
!
!! DGTSL solves a general tridiagonal linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the tridiagonal matrix.
!
!    Input/output, real ( kind = 8 ) C(N), contains the subdiagonal of the
!    tridiagonal matrix in entries C(2:N).  On output, C is destroyed.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal of the
!    matrix.  On output, D is destroyed.
!
!    Input/output, real ( kind = 8 ) E(N), contains the superdiagonal of the
!    tridiagonal matrix in entries E(1:N-1).  On output E is destroyed.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal value.
!    K, the K-th element of the diagonal becomes exactly zero.  The
!       subroutine returns if this error condition is detected.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  real ( kind = 8 ) t

  info = 0
  c(1) = d(1)

  if ( 2 <= n ) then

    d(1) = e(1)
    e(1) = 0.0D+00
    e(n) = 0.0D+00

    do k = 1, n - 1
!
!  Find the larger of the two rows.
!
      if ( abs ( c(k) ) <= abs ( c(k+1) ) ) then
!
!  Interchange rows.
!
        t = c(k+1)
        c(k+1) = c(k)
        c(k) = t

        t = d(k+1)
        d(k+1) = d(k)
        d(k) = t

        t = e(k+1)
        e(k+1) = e(k)
        e(k) = t

        t = b(k+1)
        b(k+1) = b(k)
        b(k) = t

      end if
!
!  Zero elements.
!
      if ( c(k) == 0.0D+00 ) then
        info = k
        return
      end if

      t = -c(k+1) / c(k)
      c(k+1) = d(k+1) + t * d(k)
      d(k+1) = e(k+1) + t * e(k)
      e(k+1) = 0.0D+00
      b(k+1) = b(k+1) + t * b(k)

    end do

  end if

  if ( c(n) == 0.0D+00 ) then
    info = n
    return
  end if
!
!  Back solve.
!
  b(n) = b(n) / c(n)

  if ( 1 < n ) then

    b(n-1) = ( b(n-1) - d(n-1) * b(n) ) / c(n-1)

    do k = n-2, 1, -1
      b(k) = ( b(k) - d(k) * b(k+1) - e(k) * b(k+2) ) / c(k)
    end do

  end if

  return
end