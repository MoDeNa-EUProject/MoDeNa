!> @file
!! auxilliary subroutines for calculation of state properties
!! @author    Pavel Ferkl
!! @ingroup   wall_drain
module phys_prop
    use constants, only: dp
    implicit none
    private
    integer :: Rb_kx=4,Rb_iknot=0,Rb_inbvx,Rb_nx
    real(dp), dimension(:), allocatable :: Rb_tx,Rb_coef
    integer :: por_kx=2,por_iknot=0,por_inbvx,por_nx
    real(dp), dimension(:), allocatable :: por_tx,por_coef
    integer :: visc_kx=2,visc_iknot=0,visc_inbvx,visc_nx
    real(dp), dimension(:), allocatable :: visc_tx,visc_coef
    public volume_balance,dispress,min_film_thickness,Rb_spline_ini,&
        Rb,Rb_der,porosity_spline_ini,porosity,visc_spline_ini,visc
contains
!********************************BEGINNING*************************************
!> checks whether we are losing some mass or not
pure subroutine  volume_balance(y,vt,fs)
    use globals
    real(dp), intent(out) :: vt,fs
    real(dp), dimension(:), intent(in) :: y
    integer :: i,neq
    real(dp) :: vf,vs
    neq=size(y)
    vf=0
    vs=0
    do i=1,neq
        if (y(i)<strutFilmParameter*y(1)) then
            vf=vf+2*pi*dr*(0.5_dp+i-1)*y(i)*dr
        else
            vs=vs+2*pi*dr*(0.5_dp+i-1)*y(i)*dr
        endif
    enddo
    vs=vs+pi*y(neq)/3*(rd**2+rc*rd+rc**2)-pi*y(neq)*rd**2
    vt=vf+vs
    fs=vs/vt
end subroutine volume_balance
!***********************************END****************************************


!********************************BEGINNING*************************************
!> disjoining pressure
pure subroutine dispress(h,dispr,dph)
    use globals
    real(dp), intent(in) :: h
    real(dp), intent(out) :: dispr !disjoining pressure
    real(dp), intent(out) :: dph !derivative of disjoining pressure
    dispr=bdp*((hdp/h)**(ndp-1)-(hdp/h)**(mdp-1))*(hdp-cdp)
    dph=(bdp*(hdp/h)**mdp*(hdp*mdp+cdp*(h-h*mdp))-&
        bdp*(hdp/h)**ndp*(hdp*ndp+cdp*(h-h*ndp)))/(h*hdp)
end subroutine dispress
!***********************************END****************************************


!********************************BEGINNING*************************************
!> minimum film thickness and its distance from the center
subroutine min_film_thickness(y,hmin,hloc,havg)
    use globals
    real(dp), dimension(:), intent(in) :: y
    real(dp), intent(out) :: hmin
    real(dp), intent(out) :: hloc
    real(dp), intent(out) :: havg
    integer :: i,n
    hmin=minval(y,dim=1)
    hloc=minloc(y,dim=1)*dr
    havg=0
    n=size(y,1)
    do i=1,n
        if (y(i)<strutFilmParameter*y(1)) then
            havg=havg+2*pi*dr*(0.5_dp+i-1)*y(i)*dr
        else
            exit
        endif
    enddo
    havg=havg/(pi*(dr*(0.5_dp+i-2))**2)
end subroutine min_film_thickness
!***********************************END****************************************


!********************************BEGINNING*************************************
!> time derivation of bubble radius as function of time
real(dp) function Rb_der(t)
    use bspline_module
    real(dp), intent(in) :: t
    integer :: idx,iflag
    idx=1
    call db1val(t,idx,Rb_tx,Rb_nx,Rb_kx,Rb_coef,Rb_der,iflag,Rb_inbvx)
    if (iflag /= 0) then
        print*, 'evaluation of bubble radius derivative from spline failed',&
            iflag
        stop
    endif
endfunction Rb_der
!***********************************END****************************************


!********************************BEGINNING*************************************
!> bubble radius as function of time
real(dp) function Rb(t)
    use bspline_module
    real(dp), intent(in) :: t
    integer :: idx,iflag
    idx=0
    call db1val(t,idx,Rb_tx,Rb_nx,Rb_kx,Rb_coef,Rb,iflag,Rb_inbvx)
    if (iflag /= 0) then
        print*, 'evaluation of bubble radius from spline failed',iflag
        stop
    endif
endfunction Rb
!***********************************END****************************************


!********************************BEGINNING*************************************
!> initialization of spline for bubble radius
subroutine Rb_spline_ini(bblgr_res)
    use bspline_module
    real(dp), dimension(:,:), intent(in) :: bblgr_res
    integer :: iflag
    Rb_nx=size(bblgr_res(:,1))
    if (allocated(Rb_tx)) deallocate(Rb_tx)
    if (allocated(Rb_coef)) deallocate(Rb_coef)
    allocate(Rb_tx(Rb_nx+Rb_kx),Rb_coef(Rb_nx))
    Rb_tx=0
    call db1ink(bblgr_res(:,1),Rb_nx,bblgr_res(:,2),&
        Rb_kx,Rb_iknot,Rb_tx,Rb_coef,iflag)
    Rb_inbvx=1
    if (iflag /= 0) then
        print*, 'initialization of bubble radius spline failed'
        stop
    endif
end subroutine Rb_spline_ini
!***********************************END****************************************


!********************************BEGINNING*************************************
!> porosity as function of time
real(dp) function porosity(t)
    use bspline_module
    real(dp), intent(in) :: t
    integer :: idx,iflag
    idx=0
    call db1val(t,idx,por_tx,por_nx,por_kx,por_coef,porosity,iflag,por_inbvx)
    if (iflag /= 0) then
        print*, 'evaluation of porosity from spline failed',iflag
        stop
    endif
endfunction porosity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> initialization of spline for porosity
subroutine porosity_spline_ini(bblgr_res)
    use bspline_module
    real(dp), dimension(:,:), intent(in) :: bblgr_res
    integer :: iflag
    por_nx=size(bblgr_res(:,1))
    if (allocated(por_tx)) deallocate(por_tx)
    if (allocated(por_coef)) deallocate(por_coef)
    allocate(por_tx(por_nx+por_kx),por_coef(por_nx))
    por_tx=0
    call db1ink(bblgr_res(:,1),por_nx,bblgr_res(:,3),&
        por_kx,por_iknot,por_tx,por_coef,iflag)
    por_inbvx=1
    if (iflag /= 0) then
        print*, 'initialization of porosity spline failed'
        stop
    endif
end subroutine porosity_spline_ini
!***********************************END****************************************


!********************************BEGINNING*************************************
!> viscosity as function of time
real(dp) function visc(t)
    use bspline_module
    real(dp), intent(in) :: t
    integer :: idx,iflag
    idx=0
    call db1val(t,idx,visc_tx,visc_nx,visc_kx,visc_coef,visc,iflag,visc_inbvx)
    if (iflag /= 0) then
        print*, 'evaluation of viscosity from spline failed',iflag
        stop
    endif
endfunction visc
!***********************************END****************************************


!********************************BEGINNING*************************************
!> initialization of spline for viscosity
subroutine visc_spline_ini(bblgr_res)
    use bspline_module
    real(dp), dimension(:,:), intent(in) :: bblgr_res
    integer :: iflag
    visc_nx=size(bblgr_res(:,1))
    if (allocated(visc_tx)) deallocate(visc_tx)
    if (allocated(visc_coef)) deallocate(visc_coef)
    allocate(visc_tx(visc_nx+visc_kx),visc_coef(visc_nx))
    visc_tx=0
    call db1ink(bblgr_res(:,1),visc_nx,bblgr_res(:,4),&
        visc_kx,visc_iknot,visc_tx,visc_coef,iflag)
    visc_inbvx=1
    if (iflag /= 0) then
        print*, 'initialization of viscosity spline failed'
        stop
    endif
end subroutine visc_spline_ini
!***********************************END****************************************
end module phys_prop
