!Linear interpolation from
!http://people.sc.fsu.edu/~jburkardt/f_src/pwl_interp_1d/pwl_interp_1d.html
module interpolation
    use constants, only: dp
    implicit none
    private
    public pwl_interp_1d
contains
!!********************************BEGINNING*************************************
!!minimal working example (you have to supply xy values in "refl.dat" and change value of xi)
!subroutine test_pwl_interp_1d
!    integer :: fi=154 !file index
!    integer :: i    !counter
!    integer, parameter :: nd=101   !number of data points
!    real(dp), dimension(nd) :: xd,yd   !data for interpolation
!    integer :: ni=1   !number of points, where we want to interpolate
!    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
!    real(dp) :: yi(1)   !interpolated y-values
!    open(fi,file="refl.dat")
!    do i=1,nd
!        read(fi,*) xd(i),yd(i)
!    enddo
!    close(fi)
!    xi=2e-6_dp
!    call pwl_interp_1d ( nd, xd, yd, ni, xi, yi )
!    write(*,*) xi,yi
!end subroutine test_pwl_interp_1d
!!***********************************END****************************************


subroutine pwl_interp_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************80
!
!! PWL_INTERP_1D evaluates the piecewise linear interpolant.
!
!  Discussion:
!
!    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
!    linear function which interpolates the data (XD(I),YD(I)) for I = 1
!    to ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = dp ) t
  real ( kind = dp ) xd(nd)
  real ( kind = dp ) yd(nd)
  real ( kind = dp ) xi(ni)
  real ( kind = dp ) yi(ni)

  yi(1:ni) = 0.0e+00_dp

  if ( nd == 1 ) then
    yi(1:ni) = yd(1)
    return
  end if

  do i = 1, ni

    if ( xi(i) <= xd(1) ) then

      t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
      yi(i) = ( 1.0e+00_dp - t ) * yd(1) + t * yd(2)

    else if ( xd(nd) <= xi(i) ) then

      t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
      yi(i) = ( 1.0e+00_dp - t ) * yd(nd-1) + t * yd(nd)

    else

      do k = 2, nd

        if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then

          t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
          yi(i) = ( 1.0e+00_dp - t ) * yd(k-1) + t * yd(k)
          exit

        end if

      end do

    end if

  end do

  return
end subroutine pwl_interp_1d
end module interpolation
