!> @file
!! sets up and integrates the bubble growth model
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module integration
    use constants
    use globals
    use foaming_globals_m
    use fmodena
    use modenastuff
    use in_out, only: save_integration_header,save_integration_step,&
        save_integration_close
    use model, only:odesystem,dim_var,molar_balance
    use phys_prop, only:set_initial_physical_properties
    implicit none
    private
    !time integration variables for lsode
    integer :: iout, iopt, istate, itask, itol, liw, lrw, nnz, lenrat, neq, mf
    real(dp) :: jac,tout,rtol,atol,t
    real(dp), dimension(:), allocatable :: rwork!,y
    integer, dimension(:), allocatable :: iwork
    !needed for selection of subroutine for evaluation of derivatives
    abstract interface
        subroutine sub (neq, t, y, ydot)
            use constants
            integer :: neq
            real(dp) ::  t, y(neq), ydot(neq)
        end subroutine sub
    end interface
    procedure (sub), pointer :: sub_ptr => odesystem
    public bblpreproc,bblinteg
contains
!********************************BEGINNING*************************************
!> prepares integration
subroutine bblpreproc
    write(*,*) 'preparing simulation...'
    call checks
    call set_initial_physical_properties
    call set_equation_order
    call create_mesh(0.0_dp,S0**3-R0**3,p,mshco)
    call set_initial_conditions
    call set_integrator
    write(*,*) 'done: simulation prepared'
    write(*,*)
end subroutine bblpreproc
!***********************************END****************************************


!********************************BEGINNING*************************************
!> performs some checks on the validity of input
subroutine checks
    select case (kin_model)
    case(1)
    case(3)
    case(4)
    case default
        stop 'unknown kinetic model'
    end select
end subroutine checks
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculate spatial grid points
!! point positions are stored in global array dz(p+1)
subroutine create_mesh(a,b,p,s)
    integer, intent(in) :: p ! number of inner points
    real(dp), intent(in) :: &
        a,& ! lower bound
        b,& ! upper bound
        s ! point spacing increase
	integer :: i,info
    real(dp),allocatable :: atri(:),btri(:),ctri(:),rtri(:),utri(:)
    allocate(atri(p),btri(p),ctri(p),rtri(p),utri(p),dz(p+1))
    atri=-s !lower diagonal
    btri=1+s    !main diagonal
    ctri=-1 !upper diagonal
    rtri=a  !rhs
    rtri(p)=b
    utri=rtri
    call dgtsl(p,atri,btri,ctri,utri,info)
    if ((utri(2)-utri(1))/utri(1)<10*epsilon(utri)) stop 'set smaller mesh &
        coarsening parameter'
    dz(1)=utri(1)+rtri(1)/s
    do i=2,p
        dz(i)=utri(i)-utri(i-1)
    enddo
    dz(p+1)=rtri(p)-utri(p)
    deallocate(atri,btri,ctri,rtri,utri)
end subroutine create_mesh
!***********************************END****************************************


!********************************BEGINNING*************************************
!> determine number of equations and their indexes
subroutine set_equation_order
    integer :: i
    neq=(p+1)*ngas
    if (firstrun) then
        neq = neq+4+ngas
        req=1   !radius index
        fpeq=2  !pressure index
        if (inertial_term) then
            neq=neq+1
            fpeq=fpeq+1
        endif
    else
        neq = neq+3+ngas
        fpeq=1  !pressure index
    endif
    lpeq=fpeq+ngas-1
    teq=lpeq+1   !temperature index
    xOHeq=teq+1 !polyol conversion index
    xWeq=xOHeq+1  !water conversion index
    fceq = xWeq+1    !concentration index
    if (kin_model==4) then
        allocate(kineq(20),kinsource(20))
    endif
    if (kin_model==4) then
        neq=neq+size(kineq)
        do i=1,size(kineq)
            kineq(i)=xWeq+i
        enddo
        fceq=kineq(size(kineq))+1
    endif
end subroutine set_equation_order
!***********************************END****************************************


!********************************BEGINNING*************************************
!> choose and set integrator
subroutine set_integrator
    mf=int_meth
    select case(integrator)
    case(1)
        select case(mf)
        case(10)
            allocate(rwork(20+16*neq),iwork(20))
        case(22)
            allocate(rwork(22+9*neq+neq**2),iwork(20+neq))
        case default
            stop 'unknown mf'
        end select
    case(2)
        select case(mf)
        case(10)
            allocate(rwork(20+16*neq),iwork(30))
        case(222)
            nnz=neq**2 !not sure, smaller numbers make problems for low p
            lenrat=2 !depends on dp
            allocate(rwork(int(20+(2+1._dp/lenrat)*nnz+(11+9._dp/lenrat)*neq)),&
                iwork(30))
        case default
            stop 'unknown mf'
        end select
    case default
        stop 'unknown integrator'
    end select
    itask = 1
    istate = 1
    iopt = 1
    rwork(5:10)=0
    iwork(5:10)=0
    lrw = size(rwork)
    liw = size(iwork)
    iwork(6)=maxts
    itol = 1 !don't change, or you must declare atol as atol(neq)
    rtol=rel_tol
    atol=abs_tol
end subroutine set_integrator
!***********************************END****************************************


!********************************BEGINNING*************************************
!> set initial conditions
subroutine set_initial_conditions
    integer :: i,j
    t = tstart
    tout = t+timestep
    allocate(y(neq))
    y=0
    if (firstrun) then
        y(req)=radius   !radius
        if (inertial_term) y(req+1) = 0        !velocity
    endif
    y(teq) = temp   !temperature
    y(xOHeq) = conv        !xOH
    y(xWeq) = 0        !xW
    if (kin_model==4) then
        y(kineq(1)) = 6.73000e-02_dp
        y(kineq(2)) = 1.92250e+00_dp
        y(kineq(3)) = 2.26920e+00_dp
        y(kineq(4)) = 0.00000e+00_dp
        y(kineq(5)) = 5.46200e-01_dp
        ! y(kineq(5)) = 1.0924e+00_dp
        y(kineq(6)) = 2.19790e+00_dp
        y(kineq(7)) = 1.64000e+00_dp
        y(kineq(8)) = 1.71030e+00_dp
        y(kineq(9)) = 0.00000e+00_dp
        y(kineq(10)) = 0.00000e+00_dp
        y(kineq(11)) = 0.00000e+00_dp
        y(kineq(12)) = 0.00000e+00_dp
        y(kineq(13)) = 0.00000e+00_dp
        y(kineq(14)) = 0.00000e+00_dp
        y(kineq(15)) = 0.00000e+00_dp
        y(kineq(16)) = 4.45849e+00_dp
        y(kineq(17)) = 0.00000e+00_dp
        y(kineq(18)) = 1.00000e+00_dp
        y(kineq(19)) = 60!2.27000e+01_dp
        y(kineq(20)) = 1e0_dp!8.46382e-01_dp
    endif
    do j=1,ngas
        do i=1,p+1
            y(fceq+(i-1)*ngas+j-1) = cbl(j)      !blowing agent concentration
        enddo
    enddo
    do i=1,ngas
        y(fpeq+i-1) = xgas(i+1)*(pamb+2*sigma/R0) !pressure
        if (y(fpeq+i-1)<1e-16_dp) y(fpeq+i-1)=1e-16_dp
    enddo
end subroutine set_initial_conditions
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculate growth rate
!! should be called once per timestep
subroutine growth_rate
	integer :: i
    do i=1,ngas
        grrate(i)=(mb2(i)-nold(i))/timestep
    enddo
    do i=1,ngas
        nold(i)=mb2(i)
    enddo
end subroutine growth_rate
!***********************************END****************************************


!********************************BEGINNING*************************************
!> performs integration
subroutine bblinteg
    write(*,*) 'integrating...'
    if (firstrun) then
        call save_integration_header
        call dim_var(t,y)
        call molar_balance(y)
        call growth_rate
        call save_integration_step(0)
    endif
    do iout = 1,its
        select case (integrator)
        case(1)
            call dlsode (sub_ptr, neq, y, t, tout, itol, rtol, atol, itask, &
                istate, iopt, rwork, lrw, iwork, liw, jac, mf)
        case(2)
            call dlsodes (sub_ptr, neq, y, t, tout, itol, rtol, atol, itask, &
                istate, iopt, rwork, lrw, iwork, liw, jac, mf)
        case default
            stop 'unknown integrator'
        end select
        call dim_var(t,y)
        call molar_balance(y)
        call growth_rate
        if (firstrun) call save_integration_step(iout)
        if (printout) then
            write(*,'(2x,A4,F8.3,A3,A13,F10.3,A4,A19,F8.3,A4,A9,EN12.3,A4)') &
                't = ', time, ' s,',&
                'p_b - p_o = ', bub_pres, ' Pa,', &
                'p_b - p_o - p_L = ', bub_pres-2*sigma/radius, ' Pa,',&
                'dR/dt = ', (bub_pres-2*sigma/radius)*radius/4/eta, ' m/s'
        endif
        tout = time+timestep
        if (gelpoint) exit
    end do
    if (firstrun) call save_integration_close(iout)
    write(*,*) 'done: integration'
    call destroyModenaModels
    deallocate(D,cbl,xgas,KH,Mbl,dHv,mb,mb2,mb3,avconc,pressure,&
        diff_model,sol_model,cpblg,cpbll,RWORK,IWORK,dz,y,wblpol,D0)
end subroutine bblinteg
!***********************************END****************************************
end module integration
