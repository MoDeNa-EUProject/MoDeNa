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
    use models, only: model
    use inout, only: input,output,print_header,write_header
    integer :: i, j, k, l, fi
    integer :: itol, itask, istate, iopt
    integer :: MF, ML, MU, LRW, LIW, LENRAT, NNZ, LWM, NEQ
    integer, allocatable :: IWORK(:)

    real(dp) :: tin, tout, keq
    real(dp) :: rtol, atol, jdem

    real(dp), allocatable :: ystate(:), yprime(:) ! vector of state
    real(dp), allocatable :: RWORK(:)
! -----------------------------------
! load inputs
! -----------------------------------
	call input
    do i = 1, ngas
        if ( diffModel(i) == 2  .and. modelType == "heterogeneous" ) then
            print*, "You cannot use foam diffusivity with heterogeneous model."    
            stop
        end if
    enddo
    if (sheet) then
        ncell=nint((dfoam-dsheet)/(dcell+dwall))
    else
        ncell = nint(dfoam/(dcell+dwall))
        divsheet=0
    endif
    if (modelType == "heterogeneous") then
        nfv = divsheet+ncell*(divwall+divcell)
    elseif (modelType == "homogeneous") then
        nfv = divsheet + divfoam
    endif
    neq = ngas*nfv
! -----------------------------------
! find out physical properties
! -----------------------------------
    Dgas = gasDiffusivity(temp) ! diffusivity in gas phase
    call createModels(ngas)
    eps = 1 - rhof/rhop ! porosity
    do i=1,ngas
        if (solModel(i)==1) then
            Sg(i)=Solubility(temp,i)
        endif
        if (diffModel(i)==1) then
            Dg(i)=Diffusivity(temp,i)
        endif
    enddo
    Pg = Dg*Sg*rhop/Mg/1e5_dp
    ksi = 2.0_dp
    Seff = effectiveSolubility(temp,Sg*rhop/Mg/1e5_dp,eps)
    do i = 1, ngas
        if (diffModel(i) /= 2) then
            Deff(i) = effectiveDiffusivity(dcell,dwall,Pg(i),Seff(i),ksi)
        endif
    enddo
    call print_header
    Sg=Sg*Rg*temp*rhop/(1e5*Mg)
    sheetSg=sheetSg*Rg*temp*rhop/(1e5*Mg)
! -----------------------------------
! Allocate memory for working arrays
! -----------------------------------
    allocate(ystate(1:neq))
    allocate(yprime(1:neq))
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
    if (modelType == "heterogeneous") then
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
    elseif (modelType == "homogeneous") then
        do i = 1, nfv - divsheet
            dz(k)=dfoam/(nfv - divsheet)
            mor(k)=1
            do l=1,ngas
                dif(ngas*(k-1)+l)=Deff(l)
                sol(ngas*(k-1)+l)=1
                ystate(ngas*(k-1)+l)=xg(l)*pressure*Seff(l)
            enddo
            k=k+1
        enddo
    endif
! ----------------------------------
! boundary conditions
! ----------------------------------
    if (modelType == "heterogeneous") then
        bc = pBg/Rg/temp
    elseif (modelType == "homogeneous") then
        bc = pBg*Seff
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
    itol = 1
    rtol = 1.0e-10_dp
    atol = 1.0e-6_dp
    itask = 1
    istate = 1
    iopt = 0 ! change these in case of convergence problems
! ----------------------------------
! Integration loop
! ----------------------------------
    call write_header(fi)
    call equcond(keq,ystate,ngas,nfv,mor,eps,dcell,fstrut,temp_cond)
    call output(0, 0.0_dp, ystate, neq, keq, fi)
    do i=1,nroutputs
        if (progressTime == "linear") then
            tin  = dble(i-1)*(tend-tbeg)/dble(nroutputs)+tbeg
            tout = dble(i  )*(tend-tbeg)/dble(nroutputs)+tbeg
        elseif (progressTime == "logarithmic") then
            tin  = tbeg + 10**((i - 1 - nroutputs + &
                outputsPerOrder*log10(tend - tbeg))/outputsPerOrder)
            tout = tbeg + 10**((i     - nroutputs + &
                outputsPerOrder*log10(tend - tbeg))/outputsPerOrder)
        end if
        call dlsodes(model, neq, ystate, tin, tout, itol, rtol, &
            atol, itask, istate, iopt, rwork, lrw, iwork, liw, jdem, mf)
        ! evaluating the integration
        if (istate.lt.0) then
            write(*,*) 'Something is wrong, look for ISTATE =', istate
            stop
        endif
        write(*,'(2x,A,1x,e9.3,1x,A)') 'time:', tout/(3600*24),'days'
        call equcond(keq,ystate,ngas,nfv,mor,eps,dcell,fstrut,temp_cond)
        call output(i, tout, ystate, neq, keq, fi)
    enddo
    close(fi)
    call destroyModels(ngas)
end subroutine integrate
!***********************************END****************************************
end module integration
