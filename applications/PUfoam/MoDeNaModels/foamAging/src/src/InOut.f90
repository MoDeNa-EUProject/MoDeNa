!> @file
!! subroutines for calculation of equivalent conductivity of the foam
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging

!> reads input file, packs variables to rpar and ipar variables
subroutine input(rpar, ipar)
	implicit none
!c
	integer i
	integer ipar(*)
	integer divwall, ncell
	integer nroutputs
	integer solModel(3),diffModel(3)
	double precision dcell, dwall, L
	double precision pressure ! initial conds
	double precision pBCair, pBCCO2, pBCpent 	! boundary conds
	double precision pICair, pICCO2, pICpent 	! initial conds
	double precision DO2, DN2, DCO2, Dpent, Dair, Dgas
	double precision SO2, SN2, SCO2, Spent, Sair
	double precision R, T, temp_cond, rhop
	double precision pcA, pcB, pcApcB, TcA, TcB, TcATcB
	double precision MA, MB, Mterm ,a, b, aToverTcsb
	double precision fstrut,rhof
!c
	double precision  tend,tbeg
	double precision rpar(*)					! real param

	double precision PI
	parameter (PI = 3.1415926d0)
	parameter (R = 8.314d0)

!c Read input params - pak module params :)
	open(2, file = '../foamAging.in', status = 'old')
	read(2, *) nroutputs
	read(2, *) divwall
	read(2, *) tbeg
	read(2, *) tend
	read(2, *) T
	read(2, *) temp_cond
	read(2, *) rhop
    read(2, *) pressure

	read(2, *) pBCair, pBCCO2, pBCpent  ! boundaries
	read(2, *) pICair, pICCO2, pICpent  ! initial C

	read(2, *) L
	read(2, *) dwall
	read(2, *) dcell

	read(2, *) fstrut
	read(2, *) rhof
	read(2, *) solModel
	read(2, *) Sair,SCO2,Spent
	read(2, *) diffModel
	read(2, *) Dair,DCO2,Dpent

	close(2)
	continue

	ncell = dint(L/(dcell+dwall))
    continue
	! ! computation of diffusivities and solubilities as a function of T
	! DCO2=12.3d-4*dexp(-51180.0d0/R/T)!/1e0  ! m2/s
	! DO2 =8.5d-4* dexp(-53300.0d0/R/T)
	! DN2=3.24d-3*dexp(-6927.0d0/T)
	! Dpent=1.7d-7*dexp(-4236.0d0/T)!/4e2!/2.5d0
	! Dair=(0.21d0*DO2+0.79d0*DN2)!/3e1
	!
	! ! cm3STP/cm3/Pa
	! SCO2 =  7.13d-6*T/343.0d0*dexp(-2587.0d0*(1.0d0/343.0d0-1.0d0/T))
	! Spent = 4.45d-5*T/353.0d0*dexp(-527.45d0*(1.0d0/353.0d0-1.0d0/T))
	! Sair  = -(5.0d-9 * T**2 - 4.0d-6 * T + 0.0007d0)/10
	!
	! write(*,*) 'CO2 permeability',SCO2*DCO2
	! write(*,*) 'pentane permeability',Spent*Dpent
	! write(*,*) 'air permeability',Sair*Dair

    continue
    ! gas difusivity accoriding to Bird 1975, p.505, eq. 16.3-1
    pcA = 33.5d0      !  N2
    pcB = 72.9d0      ! CO2
    pcApcB = (pcA*pcB)**(1.0d0/3.0d0) ! CO2, N2, B-1 p. 744
    TcA = 126.2d0     ! N2
    TcB = 304.2d0     ! CO2
    TcATcB = (TcA*TcB)**(5.0d0/12.0d0)
    MA = 28.02d0
    MB = 44.01d0
    Mterm = dsqrt(1/MA + 1/MB)
    a = 2.7450d-4 ! non-polar pairs
    b = 1.823d0
    aToverTcsb = a*(T/dsqrt(TcA*TcB))**b

	! pressure in atmospheres, cm2/s
    Dgas = (aToverTcsb*pcApcB*TcATcB*Mterm)*1.0d5/pressure
    continue
    Dgas = Dgas * 1.0d-4 ! m2/s
    continue

	ipar(1) = nroutputs
	ipar(2) = ncell
	ipar(4) = divwall+1		! FV in one cell
	ipar(5) = (divwall+1)*ncell	! total number of FV
	ipar(6:8) = solModel
	ipar(9:11) = diffModel

	rpar(1) = dcell
	rpar(2) = dwall
	rpar(3) = tend
	rpar(4) = Dair      ! average air
	rpar(5) = DCO2
	rpar(6) = Sair      ! average air
	rpar(7) = SCO2

	rpar(8) = pressure
	! rpar(9) = initpressure
	rpar(10)= R*T

	rpar(11) = Dpent
	rpar(12) = Spent
	rpar(13) = Dgas

    rpar(14) = pBCair
    rpar(15) = pBCCO2
    rpar(16) = pBCpent

    rpar(17) = pICair*pressure
    rpar(18) = pICCO2*pressure
    rpar(19) = pICpent*pressure

	rpar(20)=fstrut
	rpar(21)=rhof
	rpar(22)=temp_cond
	rpar(23)=rhop
	rpar(24)=tbeg

    continue

    return
end subroutine input
!c
!> saves results to file
subroutine output(iprof, time, ystate, neq)
	implicit none
	integer i, j, iprof, job
	integer nFV, onecell, ncell, neq
	integer divwall

	double precision time, test
	double precision ystate(*)
	double precision dwall, dcell, hwall

	double precision pBCair, pBCCO2, pBCpent, RT 	! boundary conds
	double precision, allocatable :: length(:)

	character(len=1) :: name_1	! one character
	character(len=2) :: name_2	! two characters
	character(len=3) :: name_3	! three characters
	character(len=4) :: name_f	! final name of file

	continue
	! ipar
    dcell = ystate(nEQ + 1)
    dwall = ystate(nEQ + 2)
    ncell   = dint(ystate(nEQ + 11))
    onecell = dint(ystate(nEQ + 12))    != dble(ipar(4))
    divwall = onecell-1
	nFV     = ncell*onecell
    hwall   = dwall/dble(divwall)

	pBCair = ystate(nEQ + 16) != rpar(14) != pBCair
    pBCCO2 = ystate(nEQ + 17) != rpar(15) != pBCCO2
    pBCpent= ystate(nEQ + 18) != rpar(16) != pBCpent
    RT     = ystate(nEQ + 10)

    if (.NOT. allocated(length)) then
        allocate (length(0:nFV))
    endif

	! compute lengths
	continue

    length(0:nFV) = 0.0d0
    do i = 0, nFV
        if (mod(i,onecell).eq.0) length(i) = length(i-1) + dcell/2.0d0
        if (mod(i,onecell).eq.1) length(i) = length(i-1) + dcell/2.0d0
        if (mod(i,onecell).gt.1) length(i) = length(i-1) + hwall
    enddo

!    continue
!     write(*,*) length(1:nFV)
     continue

    if     (iprof <  10) then
        write(name_1,'(I1)') iprof
        name_f = '000' // name_1
    elseif (iprof >= 10 .and. iprof < 100) then
        write(name_2,'(I2)') iprof
        name_f =  '00' // name_2
    elseif (iprof >= 100 .and. iprof < 1000) then
        write(name_3,'(I3)') iprof
        name_f =   '0' // name_3
    else
        write(name_f,'(I4)') iprof
    endif

    continue
	open(unit=11,file='../results/H2perm_'//trim(name_f)//'.dat')

   ! BC
	write (11,100) time/(3600.0d0*24.0d0),length(0), pBCair/RT,pBCCO2/RT,&
		pBCpent/RT
	! profiles
	do i = 1, nFV
		write (11,100) time/(3600.0d0*24.0d0),length(i),ystate(i),&
			ystate(nFV+i),ystate(2*nFV+i)
	enddo
    continue
    close (11)
	return
100   format (f8.2,F12.8,F12.6,F12.6,F12.6)
end subroutine output
!c
