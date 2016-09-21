!> @file
!! subroutines for file input and output
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module inout
	implicit none
	private
	public input,output,print_header
contains
!********************************BEGINNING*************************************
!> reads input file, packs variables to rpar and ipar variables
subroutine input()
	use constants
	use globals
	use fson
    use fson_value_m, only: fson_value_get
	type(fson_value), pointer :: json_data
	character(len=1024) :: strval
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
		json_data, "foamCondition.boundaryPressure.Cyclopentane", pBCcyp)
	call fson_get(json_data, "foamCondition.initialComposition.Air", pICair)
	call fson_get(json_data, "foamCondition.initialComposition.CO2", pICCO2)
	call fson_get(&
		json_data, "foamCondition.initialComposition.Cyclopentane", pICcyp)
	call fson_get(json_data, "morphology.foamDensity", rhof)
	call fson_get(json_data, "morphology.cellSize", dcell)
	call fson_get(json_data, "morphology.strutContent", fstrut)
	call fson_get(json_data, "morphology.wallThickness", dwall)
	call fson_get(json_data, "physicalProperties.polymerDensity", rhop)
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.Air", strval)
	if (strval=="constant") then
		solModel(1)=0
	elseif ( strval=="modena" ) then
		solModel(1)=1
	else
		print*, "Solubility model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.CO2", strval)
	if (strval=="constant") then
		solModel(2)=0
	elseif ( strval=="modena" ) then
		solModel(2)=1
	else
		print*, "Solubility model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.Cyclopentane", strval)
	if (strval=="constant") then
		solModel(3)=0
	elseif ( strval=="modena" ) then
		solModel(3)=1
	else
		print*, "Solubility model must be constant or modena"
	endif
	if (solModel(1)==0) then
		call fson_get(json_data, "physicalProperties.foam.solubility.Air", Sair)
	endif
	if (solModel(2)==0) then
		call fson_get(json_data, "physicalProperties.foam.solubility.CO2", SCO2)
	endif
	if (solModel(3)==0) then
		call fson_get(&
			json_data, "physicalProperties.foam.solubility.Cyclopentane", Scyp)
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.Air", strval)
	if (strval=="constant") then
	    diffModel(1)=0
	elseif ( strval=="modena" ) then
		diffModel(1)=1
	else
		print*, "Diffusivity model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.CO2",strval)
	if (strval=="constant") then
		diffModel(2)=0
	elseif ( strval=="modena" ) then
		diffModel(2)=1
	else
		print*, "Diffusivity model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.Cyclopentane", strval)
	if (strval=="constant") then
		diffModel(3)=0
	elseif ( strval=="modena" ) then
		diffModel(3)=1
	else
		print*, "Diffusivity model must be constant or modena"
	endif
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
			"physicalProperties.foam.diffusivity.Cyclopentane", Dcyp)
	endif
	if (sheet) then
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.Air", sheetSair)
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.CO2", sheetSCO2)
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.Cyclopentane", sheetScyp)
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.Air", sheetDair)
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.CO2", sheetDCO2)
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.Cyclopentane", sheetDcyp)
	endif
	call fson_destroy(json_data)
end subroutine input
!***********************************END****************************************


!********************************BEGINNING*************************************
!> saves results to file
subroutine output(iprof, time, ystate, neq)
	use constants
	use globals
	use model, only: ngas,dz,mor,nfv,sol
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
		write (11,100) time/(3600*24),pos,&
			ystate(ngas*(i-1)+1)*sol(ngas*(i-1)+1),&
			ystate(ngas*(i-1)+2)*sol(ngas*(i-1)+2),&
			ystate(ngas*(i-1)+3)*sol(ngas*(i-1)+3)
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


!********************************BEGINNING*************************************
subroutine print_header
	use globals
	use model, only: nfv
    print*, 'Foam:'
    write(*,'(A30,EN12.3,1x,A)') 'foam density:',rhof,'kg/m3'
    write(*,'(A30,EN12.3)') 'porosity:',eps
    write(*,'(A30,EN12.3)') 'strut content:',fstrut
    write(*,'(A30,EN12.3,1x,A)') 'cell size:',dcell,'m'
    write(*,'(A30,EN12.3,1x,A)') 'wall thickness:',dwall,'m'
    if (sheet) then
        write(*,'(A30,EN12.3,1x,A)') 'sheet thickness:',dsheet,'m'
    endif
    write(*,'(A30,EN12.3,1x,A)') 'foam thickness:',dfoam,'m'
    write(*,'(A30,I12)') 'number of cells:',ncell
    write(*,'(A30,EN12.3,1x,A)') 'aging temperature:',temp,'K'
    write(*,'(A30,EN12.3,1x,A)') 'conductivity temperature:',temp_cond,'K'
    print*, 'Physical properties:'
    write(*,'(A30,EN12.3,1x,A)') 'polymer density:',rhop,'kg/m3'
    write(*,'(A30,EN12.3,1x,A)') 'diffusivity in gas:',Dgas,'m2/s'
    write(*,'(A30,EN12.3,1x,A)') 'air diffusivity:',Dair,'m2/s'
    write(*,'(A30,EN12.3,1x,A)') 'CO2 diffusivity:',DCO2,'m2/s'
    write(*,'(A30,EN12.3,1x,A)') 'pentane diffusivity:',Dcyp,'m2/s'
	write(*,'(A30,EN12.3,1x,A)') 'air solubility:',Sair,'g/g/bar'
    write(*,'(A30,EN12.3,1x,A)') 'CO2 solubility:',SCO2,'g/g/bar'
	write(*,'(A30,EN12.3,1x,A)') 'pentane solubility:',Scyp,'g/g/bar'
	write(*,'(A30,EN12.3,1x,A)') 'air permeability:',Sair*Dair,'m2/s*g/g/bar'
    write(*,'(A30,EN12.3,1x,A)') 'CO2 permeability:',SCO2*DCO2,'m2/s*g/g/bar'
	write(*,'(A30,EN12.3,1x,A)') 'pentane permeability:',Scyp*Dcyp,&
        'm2/s*g/g/bar'
    print*, 'Numerics:'
    write(*,'(A30,I12)') 'finite volumes in wall:',divwall
    write(*,'(A30,I12)') 'finite volumes in cell:',divcell
    if (sheet) then
        write(*,'(A30,I12)') 'finite volumes in sheet:',divsheet
    endif
    write(*,'(A30,I12)') 'finite volumes in total:',nfv
    write(*,'(A30,EN12.3,1x,A)') 'initial time:',tbeg,'s'
    write(*,'(A30,EN12.3,1x,A)') 'end time:',tend,'s'
end subroutine print_header
!***********************************END****************************************
end module inout
