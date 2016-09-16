!> @file
!! Integration subroutines.
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module integration
    implicit none
    private
    public integrate
contains
!********************************BEGINNING*************************************
!> calculates the evolution of concentrations of gases
subroutine integrate
    use constants
    use globals
    use physicalProperties
    use conductivity, only: equcond
    use ioutils, only: newunit
    use model, only: modelPU,ngas,nfv,dz,dif,sol,bc,mor
    use inout, only: input,output,print_header
    integer :: i, j, k, counter, fi
	integer :: multiplicator
    integer :: itol, itask, istate, iopt
    integer :: MF, ML, MU, LRW, LIW, LENRAT, NNZ, LWM, NEQ
    integer, allocatable :: IWORK(:)

    real(dp) :: tin, tout, keq
    real(dp) :: rtol, atol, jdem

    real(dp), allocatable :: ystate(:), yprime(:) ! vector of state
    real(dp), allocatable :: RWORK(:)
    real(dp), allocatable :: yinit(:)

    ! model should be general, but physical properties and conductivity are
    ! hardcoded for ngas=3
    ngas=3
! -----------------------------------
! load inputs
! -----------------------------------
	call input
! -----------------------------------
! find out physical properties
! -----------------------------------
    Dgas= gasDiffusivity(temp)
    call createModels
    eps=1-rhof/rhop
    if (solModel(1)==1) then
        Sair=airSolubility(temp)
    endif
    if (solModel(2)==1) then
        SCO2=cdSolubility(temp)
    endif
    if (solModel(3)==1) then
        Scyp=cypSolubility(temp)
    endif
    if (diffModel(1)==1) then
        Dair=airDiffusivity(temp)
    endif
    if (diffModel(2)==1) then
        DCO2=cdDiffusivity(temp)
    endif
    if (diffModel(3)==1) then
        Dcyp=cypDiffusivity(temp)
    endif
    call print_header
    Sair=Sair*Rg*temp*1100._dp/(1e5*(0.21_dp*Mg(2)+0.79_dp*Mg(3)))
    SCO2=SCO2*Rg*temp*1100._dp/(1e5*Mg(1))
    Scyp=Scyp*Rg*temp*1100._dp/(1e5*Mg(4))
    sheetSair=sheetSair*Rg*temp*1100._dp/(1e5*(0.21_dp*Mg(2)+0.79_dp*Mg(3)))
    sheetSCO2=sheetSCO2*Rg*temp*1100._dp/(1e5*Mg(1))
    sheetScyp=sheetScyp*Rg*temp*1100._dp/(1e5*Mg(4))
! -----------------------------------
! Allocate memory for working arrays
! -----------------------------------
    if (sheet) then
        ncell=nint((dfoam-dsheet)/(dcell+dwall))
    else
        ncell = nint(dfoam/(dcell+dwall))
        divsheet=0
    endif
    nFV = divsheet+ncell*(divwall+divcell)
    nEQ = ngas*nFV
    allocate(ystate(1:nEQ))
    allocate(yprime(1:nEQ))
    allocate(dz(nfv))
    allocate(mor(nfv))
    allocate(dif(neq))
    allocate(sol(neq))
    allocate(bc(ngas))
    allocate(yinit(ngas))
! ----------------------------------
! mesh and initial conditions
! ----------------------------------
    k=1
    do i=1,divsheet
        dz(k)=dsheet/divsheet
        mor(k)=3
        dif(ngas*(k-1)+1)=sheetDair
        dif(ngas*(k-1)+2)=sheetDCO2
        dif(ngas*(k-1)+3)=sheetDcyp
        sol(ngas*(k-1)+1)=sheetSair
        sol(ngas*(k-1)+2)=sheetSCO2
        sol(ngas*(k-1)+3)=sheetScyp
        ystate(ngas*(k-1)+1)=pICair*pressure/Rg/temp
        ystate(ngas*(k-1)+2)=pICCO2*pressure/Rg/temp
        ystate(ngas*(k-1)+3)=pICcyp*pressure/Rg/temp
        k=k+1
    enddo
    do i = 1, ncell
        do j=1,divcell
            dz(k)=dcell/divcell
            mor(k)=1
            dif(ngas*(k-1)+1)=Dgas
            dif(ngas*(k-1)+2)=Dgas
            dif(ngas*(k-1)+3)=Dgas
            sol(ngas*(k-1)+1)=1
            sol(ngas*(k-1)+2)=1
            sol(ngas*(k-1)+3)=1
            ystate(ngas*(k-1)+1)=pICair*pressure/Rg/temp
            ystate(ngas*(k-1)+2)=pICCO2*pressure/Rg/temp
            ystate(ngas*(k-1)+3)=pICcyp*pressure/Rg/temp
            k=k+1
        enddo
        do j=1,divwall
            dz(k)=dwall/divwall
            mor(k)=2
            dif(ngas*(k-1)+1)=Dair
            dif(ngas*(k-1)+2)=DCO2
            dif(ngas*(k-1)+3)=Dcyp
            sol(ngas*(k-1)+1)=Sair
            sol(ngas*(k-1)+2)=SCO2
            sol(ngas*(k-1)+3)=Scyp
            ystate(ngas*(k-1)+1)=pICair*pressure/Rg/temp
            ystate(ngas*(k-1)+2)=pICCO2*pressure/Rg/temp
            ystate(ngas*(k-1)+3)=pICcyp*pressure/Rg/temp
            k=k+1
        enddo
    end do
! ----------------------------------
! boundary conditions
! ----------------------------------
    if (sheet) then
        bc(1)=pBCair/Rg/temp
        bc(2)=pBCCO2/Rg/temp
        bc(3)=pBCcyp/Rg/temp
    else
        bc(1)=pBCair/Rg/temp
        bc(2)=pBCCO2/Rg/temp
        bc(3)=pBCcyp/Rg/temp
    endif
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
    open (newunit(fi),file='../results/keq_time.out')
    write(fi,'(10A23)') '#time', 'eq.conductivity'
    call output(0, 0.0_dp, ystate, neq)
    call equcond(keq,ystate,ngas,nfv,mor,eps,dcell,fstrut,temp_cond)
    write(fi,'(10es23.15)') tbeg/(3600*24),keq*1.0e3_dp
    do i = 1, nroutputs*multiplicator       ! stabilizing multiplicator
        tin  = dble(i-1)*(tend-tbeg)/dble(nroutputs*multiplicator)+tbeg
        tout = dble(i  )*(tend-tbeg)/dble(nroutputs*multiplicator)+tbeg
100     continue    ! try to make another run for the initial step simulation
        call dlsodes(modelPU, neq, ystate, tin, tout, itol, rtol, atol, itask,&
            istate, iopt, rwork, lrw, iwork, liw, jdem, mf)
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
            call output(i/multiplicator, tout, ystate, neq)
            call equcond(keq,ystate,ngas,nfv,mor,eps,dcell,fstrut,temp_cond)
            write(fi,'(10es23.15)') tout/(3600*24),keq*1e3
        endif
    enddo
    close(10)
    close(fi)
    call destroyModels
end subroutine integrate
!***********************************END****************************************
end module integration
