!> @file
!! subroutines for growth of a single bubble and various parametric studies
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module tests
    use foaming_globals_m
    use constants
    use in_out, only:set_paths,read_inputs
    use model, only:bblpreproc,bblinteg
    implicit none
    private
    public onegrowth,eta_rm,bub_vf
contains
!********************************BEGINNING*************************************
!> simulates growth of a single bubble
subroutine onegrowth
    call set_paths
    call read_inputs
    call bblpreproc
    call bblinteg
end subroutine onegrowth
!***********************************END****************************************


!********************************BEGINNING*************************************
!> viscosity of reaction mixture as function of time
real(dp) function eta_rm(t)
    use interpolation
    real(dp) :: t
    integer :: n
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    xi(1)=t
    call pwl_interp_1d ( size(etat(:,1)), etat(:,1), etat(:,2), ni, xi, yi )
    eta_rm=yi(1)
endfunction eta_rm
!***********************************END****************************************


!********************************BEGINNING*************************************
!> volume fraction of bubbles as function of time
real(dp) function bub_vf(t)
    use interpolation
    real(dp) :: t
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    xi(1)=t
    call pwl_interp_1d ( size(port(:,1)), port(:,1), port(:,2), ni, xi, yi )
    bub_vf=yi(1)
endfunction bub_vf
!***********************************END****************************************
end module tests
