!> @file      foamConductivity/src/src/gasprop.f90
!! @ingroup   src_mod_foamConductivity
!! @author    Pavel Ferkl
!! @brief     Gas properties.
!! @details
!! Radiative properties of the gas phase.
module gasprop
    use constants
    implicit none
    private
    public abscoeffgas
contains
!********************************BEGINNING*************************************
!> Evaluates gas absorption coefficient.
real(dp) function abscoeffgas(lambda)
    use interpolation
    real(dp), intent(in) :: lambda  !< wavelength
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    xi(1)=lambda
    call pwl_interp_1d ( size(lambdagas), lambdagas, acgas, ni, xi, yi )
    abscoeffgas=yi(1)
end function abscoeffgas
!***********************************END****************************************
end module gasprop
