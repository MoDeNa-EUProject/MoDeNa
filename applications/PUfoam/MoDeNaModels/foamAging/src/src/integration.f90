!> @file      foamAging/src/src/integration.f90
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   src_mod_foamAging
!! @brief     Controls the integration.
!! @details
!! Integration is done by ODEPACK.
module integration
    implicit none
    private
    public integrate
contains
!********************************BEGINNING*************************************
!> Calculates the evolution of concentrations of gases.
!!
!! Loads the inputs, prepares the integration and holds the main integration
!! loop. And cleans up afterwards.
subroutine integrate
    use constants
    use globals
    use physicalProperties
    use conductivity, only: equcond
    use ioutils, only: newunit
    use model, only: model_heterogeneous
    use inout, only: input,output,print_header
    integer :: i, j, k, l, counter, fi
	integer :: multiplicator
    integer :: itol, itask, istate, iopt
    integer :: MF, ML, MU, LRW, LIW, LENRAT, NNZ, LWM, NEQ
    integer, allocatable :: IWORK(:)

    real(dp) :: tin, tout, keq
    real(dp) :: rtol, atol, jdem

    real(dp), allocatable :: ystate(:), yprime(:) ! vector of state
    real(dp), allocatable :: RWORK(:)
    real(dp), allocatable :: pp(:) ! partial pressure

    ! model should be general, but physical properties and conductivity are
    ! hardcoded for ngas=4
    ngas=4
    allocate(gasname(ngas))
    allocate(solModel(ngas),diffModel(ngas))
    allocate(Sg(ngas),Dg(ngas),Pg(ngas),Deff(ngas),sheetSg(ngas),sheetDg(ngas))
    allocate(pp(ngas),pBg(ngas),xg(ngas),kfoamXg(ngas),kgasXg(ngas))
    allocate(sgModena(ngas),sgInputs(ngas),sgOutputs(ngas))
    allocate(sgTemppos(ngas),sgxl1pos(ngas),sgxl2pos(ngas))
    allocate(dgModena(ngas),dgInputs(ngas),dgOutputs(ngas),dgTemppos(ngas))
    allocate(kgModena(ngas),kgInputs(ngas),kgOutputs(ngas),kgTemppos(ngas))
    ! gas names must correspond to physical properties in constants module
    gasname(1)="O2"
    gasname(2)="N2"
    gasname(3)="CO2"
    gasname(4)="CyP"
! -----------------------------------
! load inputs
! -----------------------------------
	call input
    if (sheet) then
        ncell=nint((dfoam-dsheet)/(dcell+dwall))
    else
        ncell = nint(dfoam/(dcell+dwall))
        divsheet=0
    endif
    nFV = divsheet+ncell*(divwall+divcell)
    nEQ = ngas*nFV
! -----------------------------------
! find out physical properties
! -----------------------------------
    Dgas= gasDiffusivity(temp) ! diffusivity in gas phase
    call createModels(ngas)
    eps=1-rhof/rhop
    do i=1,ngas
        if (solModel(i)==1) then
            Sg(i)=Solubility(temp,i)
        endif
        if (diffModel(i)==1) then
            Dg(i)=Diffusivity(temp,i)
        endif
    enddo
    Pg = Dg*Sg*rhop/Mg/1e5_dp
    Deff = effectiveDiffusivity(temp,dcell,dwall,Pg)
    call print_header
    Sg=Sg*Rg*temp*rhop/(1e5*Mg)
    sheetSg=sheetSg*Rg*temp*rhop/(1e5*Mg)
! -----------------------------------
! Allocate memory for working arrays
! -----------------------------------
    allocate(ystate(1:nEQ))
    allocate(yprime(1:nEQ))
    allocate(dz(nfv))
    allocate(mor(nfv))
    allocate(dif(neq))
    allocate(sol(neq))
    allocate(bc(ngas))
! ----------------------------------
! mesh and initial conditions
! ----------------------------------
    k=1
    do i=1,divsheet
        dz(k)=dsheet/divsheet
        mor(k)=3
        do l=1,ngas
            dif(ngas*(k-1)+l)=sheetDg(l)
            sol(ngas*(k-1)+l)=sheetSg(l)
            ystate(ngas*(k-1)+l)=xg(l)*pressure/Rg/temp
        enddo
        k=k+1
    enddo
    do i = 1, ncell
        do j=1,divcell
            dz(k)=dcell/divcell
            mor(k)=1
            do l=1,ngas
                dif(ngas*(k-1)+l)=Dgas
                sol(ngas*(k-1)+l)=1
                ystate(ngas*(k-1)+l)=xg(l)*pressure/Rg/temp
            enddo
            k=k+1
        enddo
        do j=1,divwall
            dz(k)=dwall/divwall
            mor(k)=2
            do l=1,ngas
                dif(ngas*(k-1)+l)=Dg(l)
                sol(ngas*(k-1)+l)=Sg(l)
                ystate(ngas*(k-1)+l)=xg(l)*pressure/Rg/temp
            enddo
            k=k+1
        enddo
    end do
! ----------------------------------
! boundary conditions
! ----------------------------------
    bc=pBg/Rg/temp
! ----------------------------------
! initialize integration
! ----------------------------------
    LENRAT = 2        ! usually for real(dp)
    NNZ = nEQ**2      ! nonzero elements  - MV CHECK
    LWM = 2*NNZ + 2*NEQ + (NNZ+10*NEQ)/LENRAT     ! MITER = 2
    LIW = 31 + NEQ + NNZ +100
    LRW = 20 + (2 + 1._dp/LENRAT)*NNZ + (11 + 9._dp/LENRAT)*NEQ
    mf = 222    ! 222 stiff, internally generated jacobi matrix
                ! 10 'normal' - chage the allocation to smaller fields
                ! 210 structure obtained NEQ+1 calls, implicit Adams,
    allocate(rwork (1:LRW))
    allocate(iwork (1:LIW))
    open(10,file='dcmp_progress.out', status='unknown')
    counter = 1
    itol = 1
    rtol = 1.0e-10_dp
    atol = 1.0e-6_dp
    itask = 1
    istate = 1
    iopt = 0
    multiplicator = 100
! ----------------------------------
! Integration loop
! ----------------------------------
    open (newunit(fi),file='keq_time.out')
    write(fi,'(10A23)') '#time', 'eq_conductivity', 'total_pressure', 'p_O2', &
        'p_N2', 'p_CO2', 'p_CP'
    call output(0, 0.0_dp, ystate, neq, pp)
    call equcond(keq,ystate,ngas,nfv,mor,eps,dcell,fstrut,temp_cond)
    write(fi,'(10es23.15)') tbeg/(3600*24),keq*1.0e3_dp,sum(pp),pp
    do i = 1, nroutputs*multiplicator       ! stabilizing multiplicator
        tin  = dble(i-1)*(tend-tbeg)/dble(nroutputs*multiplicator)+tbeg
        tout = dble(i  )*(tend-tbeg)/dble(nroutputs*multiplicator)+tbeg
100     continue    ! try to make another run for the initial step simulation
        call dlsodes(model_heterogeneous, neq, ystate, tin, tout, itol, rtol, &
            atol, itask, istate, iopt, rwork, lrw, iwork, liw, jdem, mf)
        ! evaluating the integration
        if (istate.lt.0) then
            write(*,*) 'Something is wrong, look for ISTATE =', istate
            if (istate.eq.-1) then  ! not enough steps to reach tout
                istate = 3
                iopt = 1   ! start to change something
                RWORK(5:8)=0.0e0_dp
                IWORK(5) = 0
                IWORK(6) = counter*1000
                IWORK(7) = 0
                counter = counter + 2
                write(*,*) 'MAXSTEP', IWORK(6)
                write(10,*) 'MAXSTEP', IWORK(6)
                goto 100
            endif
            ! stop
        elseif (counter>3) then
            counter=counter-1
            IWORK(6) = counter*1000
            write(*,*) 'MAXSTEP', IWORK(6)
            write(10,*) 'MAXSTEP', IWORK(6)
        endif
        ! some output
        if (mod(i,multiplicator).eq.0) then
            write(*,'(2x,A,1x,f6.1,1x,A)') 'time:', tout/(3600*24),'days'
            write(10,'(2x,A,1x,es9.3,1x,A)') 'time:', tout/(3600*24),'days'
            call output(i/multiplicator, tout, ystate, neq, pp)
            call equcond(keq,ystate,ngas,nfv,mor,eps,dcell,fstrut,temp_cond)
            write(fi,'(10es23.15)') tout/(3600*24),keq*1e3,sum(pp),pp
        endif
    enddo
    close(10)
    close(fi)
    call destroyModels(ngas)
end subroutine integrate
!***********************************END****************************************
end module integration
