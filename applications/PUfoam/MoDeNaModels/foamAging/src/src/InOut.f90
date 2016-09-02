!> @file
!! subroutines for calculation of equivalent conductivity of the foam
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module inout
	implicit none
	private
	public input,output
contains
!********************************BEGINNING*************************************
!> reads input file, packs variables to rpar and ipar variables
subroutine input(rpar, ipar)
	use constants
	use fson
    use fson_value_m, only: fson_value_get
	type(fson_value), pointer :: json_data
	integer :: i
	integer :: ipar(:)
	integer :: divwall, ncell
	integer :: nroutputs
	integer :: solModel(3),diffModel(3)
	real(dp) :: dcell, dwall, L
	real(dp) :: pressure ! initial conds
	real(dp) :: pBCair, pBCCO2, pBCpent 	! boundary conds
	real(dp) :: pICair, pICCO2, pICpent 	! initial conds
	real(dp) :: DO2, DN2, DCO2, Dpent, Dair, Dgas
	real(dp) :: SO2, SN2, SCO2, Spent, Sair
	real(dp) :: T, temp_cond, rhop
	real(dp) :: pcA, pcB, pcApcB, TcA, TcB, TcATcB
	real(dp) :: MA, MB, Mterm ,a, b, aToverTcsb
	real(dp) :: fstrut,rhof
	real(dp) ::  tend,tbeg
	real(dp) :: rpar(:)					! real param

	! Read input parameters
	json_data => fson_parse("../foamAging.json")
	call fson_get(json_data, "numberOfOutputs", nroutputs)
	call fson_get(json_data, "wallDiscretization", divwall)
	call fson_get(json_data, "timeStart", tbeg)
	call fson_get(json_data, "timeEnd", tend)
	call fson_get(json_data, "agingTemperature", T)
	call fson_get(json_data, "conductivityTemperature", temp_cond)
	call fson_get(json_data, "polymerDensity", rhop)
	call fson_get(json_data, "initialPressure", pressure)
	call fson_get(json_data, "boundaryPressure.Air", pBCair)
	call fson_get(json_data, "boundaryPressure.CO2", pBCCO2)
	call fson_get(json_data, "boundaryPressure.Cyclopentane", pBCpent)
	call fson_get(json_data, "initialComposition.Air", pICair)
	call fson_get(json_data, "initialComposition.CO2", pICCO2)
	call fson_get(json_data, "initialComposition.Cyclopentane", pICpent)
	call fson_get(json_data, "foamHalfThickness", L)
	call fson_get(json_data, "wallThickness", dwall)
	call fson_get(json_data, "cellSize", dcell)
	call fson_get(json_data, "strutContent", fstrut)
	call fson_get(json_data, "foamDensity", rhof)
	call fson_get(json_data, "solubilityModel.Air", solModel(1))
	call fson_get(json_data, "solubilityModel.CO2", solModel(2))
	call fson_get(json_data, "solubilityModel.Cyclopentane", solModel(3))
	if (solModel(1)==0) then
		call fson_get(json_data, "solubility.Air", Sair)
	endif
	if (solModel(2)==0) then
		call fson_get(json_data, "solubility.CO2", SCO2)
	endif
	if (solModel(3)==0) then
		call fson_get(json_data, "solubility.Cyclopentane", Spent)
	endif
	call fson_get(json_data, "diffusivityModel.Air", diffModel(1))
	call fson_get(json_data, "diffusivityModel.CO2", diffModel(2))
	call fson_get(json_data, "diffusivityModel.Cyclopentane", diffModel(3))
	if (diffModel(1)==0) then
		call fson_get(json_data, "diffusivity.Air", Dair)
	endif
	if (diffModel(2)==0) then
		call fson_get(json_data, "diffusivity.CO2", DCO2)
	endif
	if (diffModel(3)==0) then
		call fson_get(json_data, "diffusivity.Cyclopentane", Dpent)
	endif
	ncell = dint(L/(dcell+dwall))
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
    aToverTcsb = a*(T/dsqrt(TcA*TcB))**b

	! pressure in atmospheres, cm2/s
    Dgas = (aToverTcsb*pcApcB*TcATcB*Mterm)*1.0e5_dp/pressure
    Dgas = Dgas * 1.0e-4_dp ! m2/s

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
	rpar(10)= Rg*T

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
end subroutine input
!***********************************END****************************************


!********************************BEGINNING*************************************
!> saves results to file
subroutine output(iprof, time, ystate, neq)
	use constants
	integer :: i, j, iprof, job
	integer :: nFV, onecell, ncell, neq
	integer :: divwall

	real(dp) :: time, test
	real(dp) :: ystate(:)
	real(dp) :: dwall, dcell, hwall

	real(dp) :: pBCair, pBCCO2, pBCpent, RT 	! boundary conds
	real(dp), allocatable :: length(:)

	character(len=1) :: name_1	! one character
	character(len=2) :: name_2	! two characters
	character(len=3) :: name_3	! three characters
	character(len=4) :: name_f	! final name of file

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
    length(0:nFV) = 0.0e0_dp
    do i = 1, nFV
        if (mod(i,onecell).eq.0) length(i) = length(i-1) + dcell/2.0e0_dp
        if (mod(i,onecell).eq.1) length(i) = length(i-1) + dcell/2.0e0_dp
        if (mod(i,onecell).gt.1) length(i) = length(i-1) + hwall
    enddo
!     write(*,*) length(1:nFV)
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
	open(unit=11,file='../results/H2perm_'//trim(name_f)//'.dat')
	open(unit=12,file='../results/ppar_'//trim(name_f)//'.dat')

   ! BC
	write (11,100) time/(3600.0e0_dp*24.0e0_dp),length(0), pBCair/RT,pBCCO2/RT,&
		pBCpent/RT
	write (12,101) time/(3600.0e0_dp*24.0e0_dp),length(0), pBCair,pBCCO2,&
		pBCpent
	! profiles
	do i = 1, nFV
		write (11,100) time/(3600.0e0_dp*24.0e0_dp),length(i),ystate(i),&
			ystate(nFV+i),ystate(2*nFV+i)
	enddo
	do i = onecell, nFV, onecell
		write (12,101) time/(3600.0e0_dp*24.0e0_dp),length(i),ystate(i)*RT,&
			ystate(nFV+i)*RT,ystate(2*nFV+i)*RT
	enddo
    close(11)
	close(12)
	return
100   format (f8.2,F12.3,F12.3,F12.3,F12.3)
101   format (f8.2,F12.3,F12.3,F12.3,F12.3)
end subroutine output
!***********************************END****************************************
end module inout
