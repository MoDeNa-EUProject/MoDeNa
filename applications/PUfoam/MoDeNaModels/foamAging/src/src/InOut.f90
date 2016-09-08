!> @file
!! subroutines for file input and output
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
subroutine input()
	use constants
	use globals
	use fson
    use fson_value_m, only: fson_value_get
	type(fson_value), pointer :: json_data
	! Read input parameters
	json_data => fson_parse("../foamAging.json")
	call fson_get(json_data, "numerics.timeStart", tbeg)
	call fson_get(json_data, "numerics.timeEnd", tend)
	call fson_get(json_data, "numerics.numberOfOutputs", nroutputs)
	call fson_get(json_data, "numerics.wallDiscretization", divwall)
	call fson_get(json_data, "numerics.cellDiscretization", divcell)
	call fson_get(json_data, "numerics.sheetDiscretization", divsheet)
	call fson_get(json_data, "foamCondition.foamHalfThickness", dfoam)
	call fson_get(json_data, "foamCondition.inProtectiveSheet", sheet)
	if (sheet) then
		call fson_get(&
			json_data, "foamCondition.sheetThickness", dsheet)
	endif
	call fson_get(json_data, "foamCondition.agingTemperature", temp)
	call fson_get(json_data, "foamCondition.conductivityTemperature", temp_cond)
	call fson_get(json_data, "foamCondition.initialPressure", pressure)
	call fson_get(json_data, "foamCondition.boundaryPressure.Air", pBCair)
	call fson_get(json_data, "foamCondition.boundaryPressure.CO2", pBCCO2)
	call fson_get(&
		json_data, "foamCondition.boundaryPressure.Cyclopentane", pBCpent)
	call fson_get(json_data, "foamCondition.initialComposition.Air", pICair)
	call fson_get(json_data, "foamCondition.initialComposition.CO2", pICCO2)
	call fson_get(&
		json_data, "foamCondition.initialComposition.Cyclopentane", pICpent)
	call fson_get(json_data, "morphology.foamDensity", rhof)
	call fson_get(json_data, "morphology.cellSize", dcell)
	call fson_get(json_data, "morphology.strutContent", fstrut)
	call fson_get(json_data, "morphology.wallThickness", dwall)
	call fson_get(json_data, "physicalProperties.polymerDensity", rhop)
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.Air", solModel(1))
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.CO2", solModel(2))
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.Cyclopentane", solModel(3))
	if (solModel(1)==0) then
		call fson_get(json_data, "physicalProperties.foam.solubility.Air", Sair)
	endif
	if (solModel(2)==0) then
		call fson_get(json_data, "physicalProperties.foam.solubility.CO2", SCO2)
	endif
	if (solModel(3)==0) then
		call fson_get(&
			json_data, "physicalProperties.foam.solubility.Cyclopentane", Spent)
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.Air", diffModel(1))
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.CO2", diffModel(2))
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.Cyclopentane", diffModel(3))
	if (diffModel(1)==0) then
		call fson_get(json_data, &
			"physicalProperties.foam.diffusivity.Air", Dair)
	endif
	if (diffModel(2)==0) then
		call fson_get(json_data, &
			"physicalProperties.foam.diffusivity.CO2", DCO2)
	endif
	if (diffModel(3)==0) then
		call fson_get(json_data, &
			"physicalProperties.foam.diffusivity.Cyclopentane", Dpent)
	endif
	if (sheet) then
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.Air", sheetSair)
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.CO2", sheetSCO2)
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.Cyclopentane", sheetSpent)
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.Air", sheetDair)
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.CO2", sheetDCO2)
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.Cyclopentane", sheetDpent)
	endif
	call fson_destroy(json_data)
end subroutine input
!***********************************END****************************************


!********************************BEGINNING*************************************
!> saves results to file
subroutine output(iprof, time, ystate, neq)
	use constants
	use globals
	use model, only: ngas,dz,mor,nfv
	integer :: i, j, iprof
	integer :: neq

	real(dp) :: time,pos
	real(dp) :: ystate(:)

	character(len=1) :: name_1	! one character
	character(len=2) :: name_2	! two characters
	character(len=3) :: name_3	! three characters
	character(len=4) :: name_f	! final name of file

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

	! profiles
	do i = 1, neq/ngas
		write (11,100) time/(3600*24),pos,ystate(ngas*(i-1)+1)*Sair,&
			ystate(ngas*(i-1)+2)*SCO2,ystate(ngas*(i-1)+3)*Spent
	enddo
	do i = 1,nfv
		if (i==1) then
		    pos=dz(1)
		else
			pos=pos+(dz(i-1)+dz(i))/2
		endif
		if (mor(i)==1) then
			write (12,101) time/(3600*24),pos,ystate(ngas*(i-1)+1)*Rg*temp,&
				ystate(ngas*(i-1)+2)*Rg*temp,ystate(ngas*(i-1)+3)*Rg*temp
		endif
	enddo
    close(11)
	close(12)
100   format (f8.2,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3)
101   format (f8.2,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3)
end subroutine output
!***********************************END****************************************
end module inout
