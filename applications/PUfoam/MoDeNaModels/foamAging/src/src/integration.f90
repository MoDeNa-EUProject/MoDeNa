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
    use conductivity, only: equcond
    use physicalProperties
    use ioutils, only: newunit
    use constants
    use model, only: modelPU,initfield
    use inout, only: input,output
	implicit none
	integer ipar(20)

	integer :: nroutputs, multiplicator

	real(dp) :: rpar(24)

    integer :: itol, itask, istate, iopt
    integer :: MF, ML, MU, LRW, LIW, LENRAT, NNZ, LWM
    integer :: nFV, nEQ
    integer :: i, counter, fi
    integer, allocatable :: IWORK(:)

    real(dp) :: tin, tout, tend, keq, tbeg
    real(dp) :: rtol, atol, jdem

    real(dp), allocatable :: ystate(:), yprime(:) ! vector of state
    real(dp), allocatable :: RWORK(:)

    real(dp) :: temp,temp_cond
    real(dp) :: fstrut
    real(dp) :: rhof
    real(dp), dimension(3) :: x0
    real(dp) :: eps
    real(dp) :: rhop
    real(dp) :: ccyp

	call input(rpar, ipar)
    tbeg = rpar(24)
    temp=rpar(10)/Rg
	nroutputs = ipar(1)
	nFV = ipar(5) != (divwall+1)*ncell	! total number of FV
    nEQ = 3*nFV
    solModel=ipar(6:8)
    diffModel = ipar(9:11)
! -----------------------------------
! find out physical properties
! -----------------------------------
    call createModels
    fstrut=rpar(20)
    rhof=rpar(21)
    rhop=rpar(23)
    eps=1-rhof/rhop
    temp_cond=rpar(22)
    if (ipar(6)==1) then
        rpar(6)=airSolubility(temp)
    endif
    if (ipar(7)==1) then
        rpar(7)=cdSolubility(temp)
    endif
    if (ipar(8)==1) then
        rpar(12)=cypSolubility(temp)
    endif
    if (ipar(9)==1) then
        rpar(4)=airDiffusivity(temp)
    endif
    if (ipar(10)==1) then
        rpar(5)=cdDiffusivity(temp)
    endif
    if (ipar(11)==1) then
        rpar(11)=cypDiffusivity(temp)
    endif
	write(*,*) 'air solubility',rpar(6)
    write(*,*) 'CO2 solubility',rpar(7)
	write(*,*) 'pentane solubility',rpar(12)
	write(*,*) 'air diffusivity',rpar(4)
    write(*,*) 'CO2 diffusivity',rpar(5)
	write(*,*) 'pentane diffusivity',rpar(11)
	write(*,*) 'air permeability',rpar(6)*rpar(4)
    write(*,*) 'CO2 permeability',rpar(7)*rpar(5)
	write(*,*) 'pentane permeability',rpar(12)*rpar(11)
    ! stop
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
    allocate( ystate(1:nEQ  + 20)   )
    allocate( yprime(1:nEQ)   )
    allocate( rwork (1:LRW) )
    allocate( iwork (1:LIW) )
!c ----------------------------------
!c pack ystate
!c ----------------------------------
	ystate(nEQ + 1 ) = rpar(1) != dcell
	ystate(nEQ + 2 ) = rpar(2) != dwall
	ystate(nEQ + 3 ) = rpar(3) != tend
	ystate(nEQ + 4 ) = rpar(4) != 0.21e0_dp*DO2+0.79e0_dp*DN2      ! average air
	ystate(nEQ + 5 ) = rpar(5) != DCO2
	ystate(nEQ + 6 ) = rpar(6) != 0.21e0_dp*SO2+0.79e0_dp*SN2      ! average air
	ystate(nEQ + 7 ) = rpar(7) != SCO2
	ystate(nEQ + 8 ) = rpar(8) != pressure
	! ystate(nEQ + 9 ) = rpar(9) != initpressure
	ystate(nEQ + 10) = rpar(10)!= R*T
	ystate(nEQ + 11) = dble(ipar(2))! ncell
	ystate(nEQ + 12) = dble(ipar(4))  ! onecell
    ystate(nEQ + 13) = rpar(11) != Dpent
	ystate(nEQ + 14) = rpar(12) != Spent
	ystate(nEQ + 15) = rpar(13) != Dgas

	ystate(nEQ + 16) = rpar(14) != pBCair
    ystate(nEQ + 17) = rpar(15) != pBCCO2
    ystate(nEQ + 18) = rpar(16) != pBCpent
!c ----------------------------------
! initial conditions
!c ----------------------------------
    call initfield(ystate, nFV, nEQ, rpar)
!c ----------------------------------
! initialize integration
!c ----------------------------------
    open(10,file='dcmp_progress.out', status='unknown')
    tend = rpar(3)
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
