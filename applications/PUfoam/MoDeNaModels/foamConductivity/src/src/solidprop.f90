!> @file
!! subroutines for evaluation of radiative properties of solid
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module solidprop
    use constants
    implicit none
    private
    public optconst
contains
!********************************BEGINNING*************************************
!> determine optical constants n,k for one wavelength
subroutine optconst(lambda,n,k)
    use quadpack
    use interpolation
    real(dp), intent(in) :: lambda  !wavelength
    real(dp), intent(out) :: n,k  !complex refractive index
    !qags variables
    real(dp) :: a !start point of integration
    real(dp) :: abserr
    real(dp) :: b !end point of integration
    real(dp), parameter :: epsabs = 0.0e0_dp
    real(dp), parameter :: epsrel = 0.001e0_dp
    integer :: ier
    integer :: neval
    real(dp) :: res
    !end of qags variables
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    a=lambdan(1)
    b=lambdan(size(lambdan))
    if (lambda<a) then
        write(*,*) 'No data for such low wavelength.'
        write(*,*) 'Minimum wavelength is',a
        stop
    elseif (lambda>b) then
        call qag ( nwew, a, b, epsabs, epsrel, 1, res, abserr, neval, ier )
        if (ier /= 0) then
            write(*,*) 'qag returned',ier
            write(*,*) 'wavelength',lambda
            write(*,*) 'new real part of refractive index not calculated'
            stop
        endif
        n=res
        call qag ( kwew, a, b, epsabs, epsrel, 1, res, abserr, neval, ier )
        if (ier /= 0) then
            write(*,*) 'qag returned',ier
            write(*,*) 'wavelength',lambda
            write(*,*) 'new imaginery part of refractive index not calculated'
            stop
        endif
        k=res
        call qag ( Planck2, a, b, epsabs, epsrel, 1, res, abserr, neval, ier )
        if (ier /= 0) then
            write(*,*) 'qag returned',ier
            write(*,*) 'wavelength',lambda
            write(*,*) 'new refractive index not calculated'
            stop
        endif
        n=n/res
        k=k/res
    else
        xi(1)=lambda
        call pwl_interp_1d ( size(lambdan), lambdan, nwl, ni, xi, yi )
        n=yi(1)
        call pwl_interp_1d ( size(lambdak), lambdak, kwl, ni, xi, yi )
        k=yi(1)
    endif
end subroutine optconst
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates real part of index of refraction times emissive power
real(dp) function nwew ( lambda )
!***************************DECLARATION******************************
    use interpolation
    real(dp) :: lambda
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
!******************************BODY**********************************
    xi(1)=lambda
    call pwl_interp_1d ( size(lambdan), lambdan, nwl, ni, xi, yi )
    nwew=yi(1)*Planck(tmean,lambda)
end function nwew
!***********************************END****************************************


!********************************BEGINNING*************************************
!> evaluates real part of index of refraction times emissive power
real(dp) function kwew ( lambda )
!***************************DECLARATION******************************
    use interpolation
    real(dp) :: lambda
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
!******************************BODY**********************************
    xi(1)=lambda
    call pwl_interp_1d ( size(lambdak), lambdak, kwl, ni, xi, yi )
    kwew=yi(1)*Planck(tmean,lambda)
end function kwew
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determines emissive power
real(dp) function Planck(temp,lambda)
!***************************DECLARATION******************************
    real(dp), intent(in) :: temp    !temperature
    real(dp), intent(in) :: lambda  !wavelength
    real(dp) :: n
!******************************BODY**********************************
    n=1.57_dp !just guess, equation for emissive power is correct only
    ! for constant n
    Planck=2*pi*hPc*c0**2/(n**2*lambda**5*(exp(hPc*c0/(n*lambda*kb*temp))-1))
end function Planck
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determines emissive power
real(dp) function Planck2(lambda)
!***************************DECLARATION******************************
    real(dp), intent(in) :: lambda  !wavelength
!******************************BODY**********************************
    Planck2=Planck(tmean,lambda)
end function Planck2
!***********************************END****************************************
end module solidprop
