!> @file
!! subroutines for growth of a single bubble and various parametric studies
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module tests
    use foaming_globals_m
    use constants
    use in_out, only:set_paths,read_inputs,load_old_results,&
        read_inputs_json
    use integration, only:bblpreproc,bblinteg
    implicit none
    private
    public onegrowth,eta_rm,bub_vf,secondgrowth
contains
!********************************BEGINNING*************************************
!> simulates growth of a single bubble
!! you must set firstrun variable before calling this subroutine
!! if firstrun==.true.:
!!     no need to set anything
!! if firstrun==.false.:
!!     set tend
!!     set bub_rad
!!     set bub_inx
subroutine onegrowth
    call set_paths
    call read_inputs
    call bblpreproc
    call bblinteg
end subroutine onegrowth
!***********************************END****************************************


!********************************BEGINNING*************************************
!> simulates growth of a single bubble
!! uses precalculated evolution of bubble radius to calculate bubble pressure
subroutine secondgrowth
    use globals
    ! prepare to work with files
    call set_paths
    ! obtain tend and clean after yourself
    firstrun=.true.
    call read_inputs_json
    deallocate(D,cbl,xgas,KH,Mbl,&
        dHv,mb,mb2,mb3,avconc,pressure,&
        diff_model,sol_model,cpblg,cpbll,&
        wblpol,D0)
    ! feed results to function Rb(t)
    call load_old_results
    bub_inx=1
    ! simulation
    firstrun=.false.
    call onegrowth
end subroutine secondgrowth
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
