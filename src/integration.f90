!sets up and integrates the models
module integration
    use globals
    implicit none
    private
    public preprocess,integrate
    !time integration variables for lsode
    real(dp), dimension(:), allocatable :: rwork,y
    integer, dimension(:), allocatable :: iwork
    real(dp) :: jac,tout,rtol,atol,t
    integer :: iopt, istate, itask, itol, liw, lrw, nnz, lenrat, mf, neq
    !needed for selection of model subroutine
    abstract interface
        subroutine sub (neq, t, y, ydot)
            use constants
            integer :: neq
            real(dp) ::  t, y(neq), ydot(neq)
        end subroutine sub
    end interface
    procedure (sub), pointer :: odesystem_ptr => null()
    real(dp) :: vt0,told
    real(dp), dimension(:), allocatable :: yold
contains
!********************************BEGINNING*************************************
! prepares integration
subroutine preprocess
    use in_out, only: read_inputs
    use phys_prop, only: volume_balance
    use model, only: odesystem,set_initial_conditions
    real(dp) :: vf,vs,vt
    write(*,*) 'wellcome to wall drainage.'
    call read_inputs
    odesystem_ptr=>odesystem
    neq=meshpoints
    mf=int_method
    t=initialTime
    atol=int_abstol
    rtol=int_reltol
    q=1.0e-15_dp
    s=1/sqrt(3._dp)
    allocate(y(neq),yold(neq))
    nnz=neq**2 !i really don't know, smaller numbers can make problems
    lenrat=2 !depends on dp
    allocate(rwork(int(20+(2+1._dp/lenrat)*nnz+(11+9._dp/lenrat)*neq)),&
        iwork(30))
    itask = 1
    istate = 1
    iopt = 1
    rwork(5:10)=0
    iwork(5:10)=0
    lrw = size(rwork)
    liw = size(iwork)
    iwork(6)=maxts
    tout =t+timestep
    itol = 1 !don't change, or you must declare atol as atol(neq)
    call set_initial_conditions(y)
    call volume_balance(y,vf,vs,vt)
    vt0=vt
end subroutine preprocess
!***********************************END****************************************


!********************************BEGINNING*************************************
! prepares integration
subroutine integrate
    use Solve_NonLin, only: hbrd
    use in_out, only: save_int_header,save_int_step,save_int_close
    integer :: i,fi,fi2
    integer, parameter :: n=1
    integer :: info
    real (dp), dimension(n) :: x,fvec,diag
    call save_int_header(y,t,dr)
    do i=1,its
        told=t
        yold=y
        x(1)=q
        ! at each time step find the flux, which will preserve total volume
        ! of the system (film + strut)
        call hbrd(drain_residual,n,x,fvec,epsilon(pi),ae_tol,info,diag)
        if (info /= 1) then
            print*, 'Flux not found.'
            print*, 'Hbrd returned info = ',info
            stop
        endif
        call save_int_step(y,t,dr)
        tout = tout+timestep
    enddo
    call save_int_close
    write(*,*) 'program exited normally.'
end subroutine integrate
!***********************************END****************************************


!********************************BEGINNING*************************************
!> residual function for the draininng
subroutine drain_residual(n,x,fvec,iflag)
    use phys_prop, only: volume_balance
    integer, intent(in) :: n
    integer, intent(inout) :: iflag
    real(dp), dimension(n), intent(in) :: x
    real(dp), dimension(n), intent(out) :: fvec
    real(dp) :: vf,vs,vt
    q=x(1)
    t=told
    y=yold
    istate=1
    call dlsodes (odesystem_ptr, neq, y, t, tout, itol, rtol, atol, itask, &
        istate, iopt, rwork, lrw, iwork, liw, jac, mf)
    rc=rc0+gr*t
    rd=rc-y(neq)/sqrt(3.0_dp)
    dr=rd/neq
    call volume_balance(y,vf,vs,vt)
    fvec(1)=sqrt((vt-vt0)**2)
    print*, x(1),fvec(1)
end subroutine drain_residual
!***********************************END****************************************
end module integration
