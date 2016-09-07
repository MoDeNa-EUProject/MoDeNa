!> @file
!! Integration subroutines.
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module integration
    implicit none
    private
    public degas
contains
!********************************BEGINNING*************************************
!> calculates the evolution of concentrations of gases
subroutine degas
    use constants
    use globals
    use physicalProperties
    use conductivity, only: equcond
    use ioutils, only: newunit
    use model, only: modelPU,initfield,ngas,nfv,dz,dif,sol,bc
    use inout, only: input,output
	integer :: multiplicator
    integer :: itol, itask, istate, iopt
    integer :: MF, ML, MU, LRW, LIW, LENRAT, NNZ, LWM, NEQ
    integer :: i, j, k, counter, fi
    integer, allocatable :: IWORK(:)

    real(dp) :: tin, tout, keq
    real(dp) :: rtol, atol, jdem

    real(dp), allocatable :: ystate(:), yprime(:) ! vector of state
    real(dp), allocatable :: RWORK(:)
    real(dp), allocatable :: yinit(:)

    real(dp) :: eps
	real(dp) :: pcA, pcB, pcApcB, TcA, TcB, TcATcB
	real(dp) :: MA, MB, Mterm ,a, b, aToverTcsb

    ! model should be general, but physical properties hardcoded for ngas=3
    ngas=3
    ! load inputs
	call input
    if (sheet) then
        ncell=nint((dfoam-dsheet)/(dcell+dwall))
    else
        ncell = nint(dfoam/(dcell+dwall))
        divsheet=0
    endif
    ! gas difusivity accoriding to Bird 1975, p.505, eq. 16.3-1
    pcA = 33.5e0_dp      !  N2
    pcB = 72.9e0_dp      ! CO2
    pcApcB = (pcA*pcB)**(1.0e0_dp/3.0e0_dp) ! CO2, N2, B-1 p. 744
    TcA = 126.2e0_dp     ! N2
    TcB = 304.2e0_dp     ! CO2
    TcATcB = (TcA*TcB)**(5.0e0_dp/12.0e0_dp)
    MA = 28.02e0_dp
    MB = 44.01e0_dp
    Mterm = dsqrt(1/MA + 1/MB)
    a = 2.7450e-4_dp ! non-polar pairs
    b = 1.823e0_dp
    aToverTcsb = a*(temp/dsqrt(TcA*TcB))**b
    ! pressure in atmospheres, cm2/s
    Dgas = (aToverTcsb*pcApcB*TcATcB*Mterm)*1.0e5_dp/pressure
    Dgas = Dgas * 1.0e-4_dp ! m2/s
	nFV = divsheet+ncell*(divwall+divcell)
    nEQ = ngas*nFV
! -----------------------------------
! find out physical properties
! -----------------------------------
    call createModels
    eps=1-rhof/rhop
    if (solModel(1)==1) then
        Sair=airSolubility(temp)
    endif
    if (solModel(2)==1) then
        SCO2=cdSolubility(temp)
    endif
    if (solModel(3)==1) then
        Spent=cypSolubility(temp)
    endif
    if (diffModel(1)==1) then
        Dair=airDiffusivity(temp)
    endif
    if (diffModel(2)==1) then
        DCO2=cdDiffusivity(temp)
    endif
    if (diffModel(3)==1) then
        Dpent=cypDiffusivity(temp)
    endif
	write(*,*) 'air solubility',Sair
    write(*,*) 'CO2 solubility',SCO2
	write(*,*) 'pentane solubility',Spent
	write(*,*) 'air diffusivity',Dair
    write(*,*) 'CO2 diffusivity',DCO2
	write(*,*) 'pentane diffusivity',Dpent
	write(*,*) 'air permeability',Sair*Dair
    write(*,*) 'CO2 permeability',SCO2*DCO2
	write(*,*) 'pentane permeability',Spent*Dpent
!c -----------------------------------
!c Allocate memory for working arrays
!c -----------------------------------
    LENRAT = 2        ! usually for real(dp)
    NNZ = nEQ**2      ! nonzero elements  - MV CHECK
    LWM = 2*NNZ + 2*NEQ + (NNZ+10*NEQ)/LENRAT     ! MITER = 2
    LIW = 31 + NEQ + NNZ +100
    LRW = 20 + (2 + 1._dp/LENRAT)*NNZ + (11 + 9._dp/LENRAT)*NEQ
    mf = 222    ! 222 stiff, internally generated jacobi matrix
                ! 10 'normal' - chage the allocation to smaller fields
                ! 210 structure obtained NEQ+1 calls, implicit Adams,
!c -----------------------------------
!c Allocate memory for working arrays
!c -----------------------------------
    allocate(ystate(1:nEQ))
    allocate(yprime(1:nEQ))
    allocate(rwork (1:LRW))
    allocate(iwork (1:LIW))
    allocate(dz(nfv))
    allocate(dif(neq))
    allocate(sol(neq))
    allocate(bc(ngas))
    allocate(yinit(ngas))
!c ----------------------------------
!c mesh
!c ----------------------------------
    k=1
    do i=1,divsheet
        dz(k)=dsheet/divsheet
        dif(ngas*(k-1)+1)=sheetDair
        dif(ngas*(k-1)+2)=sheetDCO2
        dif(ngas*(k-1)+3)=sheetDpent
        sol(ngas*(k-1)+1)=sheetSair
        sol(ngas*(k-1)+2)=sheetSCO2
        sol(ngas*(k-1)+3)=sheetSpent
        k=k+1
    enddo
    do i = 1, ncell
        do j=1,divwall
            dz(k)=dwall/divwall
            dif(ngas*(k-1)+1)=Dair
            dif(ngas*(k-1)+2)=DCO2
            dif(ngas*(k-1)+3)=Dpent
            sol(ngas*(k-1)+1)=Sair
            sol(ngas*(k-1)+2)=SCO2
            sol(ngas*(k-1)+3)=Spent
            k=k+1
        enddo
        do j=1,divcell
            dz(k)=dcell/divcell
            dif(ngas*(k-1)+1)=Dgas
            dif(ngas*(k-1)+2)=Dgas
            dif(ngas*(k-1)+3)=Dgas
            sol(ngas*(k-1)+1)=1
            sol(ngas*(k-1)+2)=1
            sol(ngas*(k-1)+3)=1
            k=k+1
        enddo
    end do
!c ----------------------------------
! initial conditions
!c ----------------------------------
    yinit(1)=pICair*pressure/Sair
    yinit(2)=pICCO2*pressure/SCO2
    yinit(3)=pICpent*pressure/Spent
    call initfield(neq,ystate,yinit) !TODO: wrong in gas and sheet
    print*, ystate
    stop
!c ----------------------------------
! boundary conditions
!c ----------------------------------
    bc(1)=pBCair/Sair
    bc(2)=pBCCO2/SCO2
    bc(3)=pBCpent/Spent
!c ----------------------------------
! initialize integration
!c ----------------------------------
    open(10,file='dcmp_progress.out', status='unknown')
    counter = 1
    itol = 1
    rtol = 0
    atol = 1.0e-6_dp

    itask = 1
    istate = 1
    iopt = 0

    multiplicator = 100
!c ----------------------------------
!c Integration loop
!c ----------------------------------
    open (newunit(fi),file='../results/keq_time.out')
    write(fi,'(10A23)') '#time', 'eq.conductivity'
    call output(0, 0.0_dp, ystate, neq)
    call equcond(keq,ystate,neq,eps,fstrut,temp_cond)
    write(fi,'(10es23.15)') tbeg/(3600.0e0_dp*24.0e0_dp),keq*1.0e3_dp
    do i = 1, nroutputs*multiplicator       ! stabilizing multiplicator
        tin = dble(i-1)*(tend -tbeg)/dble(nroutputs*multiplicator)+tbeg
        tout = dble(i  )*(tend-tbeg)/dble(nroutputs*multiplicator)+tbeg
100     continue    ! try to make another run for the initial step simulation
        call dlsodes (modelPU, neq, ystate, tin, tout, itol, rtol, atol, itask,&
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
        elseif (counter>3) then
            counter=counter-1
            IWORK(6) = counter*1000
            write(*,*) 'MAXSTEP', IWORK(6)
            write(10,*) 'MAXSTEP', IWORK(6)
        endif
        ! some output
        if (mod(i,multiplicator).eq.0) then
            write(*,'(2x,A,1x,f6.1,1x,A)') &
                'time:', tout/(3600.0e0_dp*24.0e0_dp),'days'
            write(10,'(2x,A,1x,es9.3,1x,A)') &
                'time:', tout/(3600.0e0_dp*24.0e0_dp),'days'
            call output(i/multiplicator, tout, ystate, neq)
            call equcond(keq,ystate,neq,eps,fstrut,temp_cond)
            write(fi,'(10es23.15)') tout/(3600.0e0_dp*24.0e0_dp),keq*1e3
        endif
    enddo
    close(10)
    close(fi)
    call destroyModels
end subroutine degas
!***********************************END****************************************
end module integration
