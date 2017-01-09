!> @file      bubbleGrowth/src/src/tests.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_bubbleGrowth
!! @brief     Top level subroutines.
!! @details
!! Subroutines for growth of a single bubble and various parametric studies.
!! Also subroutines, which we want to make public in case the model is compiled
!! into a library.
module tests
    use foaming_globals_m
    use constants
    use in_out, only:set_paths,read_inputs,load_old_results
    use integration, only:bblpreproc,bblinteg
    implicit none
    private
    public onegrowth,eta_rm,bub_vf,secondgrowth,shooting_method,&
        shooting_method_test
contains
!********************************BEGINNING*************************************
!> Simulates growth of a single bubble.
!!
!! You must set firstrun variable before calling this subroutine.\n
!! firstrun==.true. => no need to set anything else\n
!! firstrun==.false.: => set tend, bub_ra, bub_inx
subroutine onegrowth
    use phys_prop, only:Rb_initialized
    Rb_initialized=.false.
    call set_paths
    call read_inputs
    call bblpreproc
    call bblinteg
    Rb_initialized=.false.
end subroutine onegrowth
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Simulates growth of a single bubble.
!!
!! Uses precalculated evolution of bubble radius to calculate bubble pressure.
!! Does not use momentum balance in the system of equations.
subroutine secondgrowth
    use globals
    use fson
    use fson_value_m, only: fson_value_get
    real(dp) :: tt,rr
    type(fson_value), pointer :: json_data
    ! prepare to work with files
    call set_paths
    ! obtain tend
    write(*,*) 'loading input file ',TRIM(inputs)
    json_data => fson_parse(inputs)
    call fson_get(json_data, "bubbleGrowth.finalTime", tend)
    ! feed results to function Rb(t)
    call load_old_results
    bub_inx=1
    ! simulation
    firstrun=.false.
    call onegrowth
end subroutine secondgrowth
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Function calculates bubble radius based on initial bubble radius.
!!
!! Sets the initial radius to first element of \p rad_ini, makes the simulation
!! and returns the bubble radius at the end.
real(dp) function nextRadius(n,rad_ini)
    use globals
    integer, intent(in) :: n !< length of \p rad_ini
    real(dp), dimension(n), intent(in) :: rad_ini !< array holding the
    !! previously calculated evolution of bubble radius
    R0=rad_ini(1)
    call onegrowth
    nextRadius=radius
end function nextRadius
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Residual function for the shooting method.
!!
!! Residual function for the Powell method as implemented in
!! @ref bubbleGrowth/src/src/hbrd.f90
subroutine rad_residual(n,x,fvec,iflag)
    use globals
    integer, intent(in) :: n !< number of equations
    integer, intent(inout) :: iflag !< should not be changed in this function
    real(dp), dimension(n), intent(in) :: x !< array with independent variables
    real(dp), dimension(n), intent(out) :: fvec !< function values
    fvec(1)=sqrt((nextRadius(n,x)-goalRadius)**2)
    print*, x(1),fvec(1)
end subroutine rad_residual
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Finds foaming conditions that end with desired bubble radius.
!!
!! Uses shooting method to find the initial bubble radius, which leads to
!! the desired final radius.
subroutine shooting_method
    use globals
    use Solve_NonLin
    use fson
    use fson_value_m, only: fson_value_get
    integer, parameter :: n=1
    integer :: info
    real (dp) :: rad_ini,frad_ini,tol=1e-8_dp
    real (dp), dimension(n) :: x,fvec,diag
    type(fson_value), pointer :: json_data
    ! prepare to work with files
    call set_paths
    ! obtain tend
    write(*,*) 'loading input file ',TRIM(inputs)
    json_data => fson_parse(inputs)
    call fson_get(json_data, "initialConditions.bubbleRadius", R0)
    rad_ini=R0
    goalRadius=bub_rad(size(bub_rad(:,1)),bub_inx+1)
    firstrun=.true.
    shooting=.true.
    x(1)=rad_ini
    call hbrd(rad_residual,n,x,fvec,epsilon(pi),tol,info,diag)
    if (info /= 1) then
        write(unit=*, fmt=*) 'Shooting method failed.'
        write(unit=*, fmt=*) 'Hbrd returned info = ',info
        stop
    endif
    R0=x(1)
    call onegrowth
    write(unit=*, fmt=*) 'First guess initial radius = ',rad_ini
    write(unit=*, fmt=*) 'Converged initial radius = ',R0
    write(unit=*, fmt=*) 'Goal radius = ',goalRadius
    write(unit=*, fmt=*) 'Final converged radius = ',radius
end subroutine shooting_method
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Test for the shooting_method.
!!
!! Runs the simulation and then tries to determine initial radius, which would
!! lead to a slightly larger final radius.
subroutine shooting_method_test
    use globals
    use Solve_NonLin
    integer, parameter :: n=1
    integer :: info
    real (dp) :: rad_ini,frad_ini,tol=1e-8_dp
    real (dp), dimension(n) :: x,fvec,diag
    ! first simulation
    firstrun=.true.
    shooting=.false.
    call onegrowth
    rad_ini=R0
    frad_ini=radius
    goalRadius=radius*(1+1.e-4_dp)
    shooting=.true.
    x(1)=rad_ini
    call hbrd(rad_residual,n,x,fvec,epsilon(pi),tol,info,diag)
    if (info /= 1) then
        write(unit=*, fmt=*) 'Shooting method failed.'
        write(unit=*, fmt=*) 'Hbrd returned info = ',info
        stop
    endif
    R0=x(1)
    call onegrowth
    write(unit=*, fmt=*) 'First guess initial radius = ',rad_ini
    write(unit=*, fmt=*) 'Converged initial radius = ',R0
    write(unit=*, fmt=*) 'Goal radius = ',goalRadius
    write(unit=*, fmt=*) 'First guess final radius = ',frad_ini
    write(unit=*, fmt=*) 'Final converged radius = ',radius
end subroutine shooting_method_test
!***********************************END****************************************



!********************************BEGINNING*************************************
!> Viscosity of reaction mixture as function of time.
!!
!! Uses linear interpolation on \p etat array.
real(dp) function eta_rm(t)
    use interpolation
    real(dp), intent(in) :: t !< time
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
!> Volume fraction of bubbles as function of time.
!!
!! Uses linear interpolation on \p port array.
real(dp) function bub_vf(t)
    use interpolation
    real(dp), intent(in) :: t !< time
    integer :: ni=1   !number of points, where we want to interpolate
    real(dp) :: xi(1)   !x-values of points, where we want to interpolate
    real(dp) :: yi(1)   !interpolated y-values
    xi(1)=t
    call pwl_interp_1d ( size(port(:,1)), port(:,1), port(:,2), ni, xi, yi )
    bub_vf=yi(1)
endfunction bub_vf
!***********************************END****************************************
end module tests
