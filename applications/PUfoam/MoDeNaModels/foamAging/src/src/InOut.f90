!> @file      foamAging/src/src/InOut.f90
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   src_mod_foamAging
!! @brief     File input and output.
!! @details
!! Inputs are loaded from JSON file, pressure profiles are saved in text files.
module inout
	implicit none
	private
	public input,output,print_header
contains
!********************************BEGINNING*************************************
!> Reads input file.
!!
!! Inputs are saved to global variables.
subroutine input()
	use constants
	use globals
	use physicalProperties
	use ioutils
	use fson
    use fson_value_m, only: fson_value_get
	type(fson_value), pointer :: json_data
	character(len=1024) :: strval
    character(len=99) :: after_foaming,after_foaming0='after_foaming.txt'
    character(len=99) :: bg_res='../../foamExpansion/results/bubbleGrowth/'
    character(len=99) :: qmom0D_res='../../foamExpansion/results/CFD0D/'
    character(len=99) :: qmom3D_res='../../foamExpansion/results/CFD3D/'
	real(dp) :: matr(7)
	integer :: fi
	! Read input parameters
	json_data => fson_parse("../inputs/foamAging.json")
	call fson_get(json_data, "modelType", modelType)
	call fson_get(json_data, "foamCondition.inProtectiveSheet", sheet)
	call fson_get(json_data, "numerics.timeStart", tbeg)
	call fson_get(json_data, "numerics.timeEnd", tend)
	call fson_get(json_data, "numerics.numberOfOutputs", nroutputs)
	call fson_get(json_data, "numerics.wallDiscretization", divwall)
	call fson_get(json_data, "numerics.cellDiscretization", divcell)
	if (sheet) then
		call fson_get(json_data, "numerics.sheetDiscretization", divsheet)
	endif
	call fson_get(json_data, "foamCondition.foamHalfThickness", dfoam)
	if (sheet) then
		call fson_get(&
			json_data, "foamCondition.sheetThickness", dsheet)
	endif
	call fson_get(json_data, "foamCondition.agingTemperature", temp)
	call fson_get(json_data, "foamCondition.conductivityTemperature", temp_cond)
	call fson_get(json_data, "foamCondition.initialPressure", pressure)
	call fson_get(json_data, "foamCondition.boundaryPressure.O2", pBg(1))
	call fson_get(json_data, "foamCondition.boundaryPressure.N2", pBg(2))
	call fson_get(json_data, "foamCondition.boundaryPressure.CO2", pBg(3))
	call fson_get(&
		json_data, "foamCondition.boundaryPressure.Cyclopentane", pBg(4))
	call fson_get(json_data, "sourceOfProperty.gasComposition", strval)
    if (strval=="BubbleGrowth") then
        after_foaming=TRIM(ADJUSTL(bg_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xg(1)=0 ! no oxygen
		xg(2)=0 ! no nitrogen
        xg(3)=matr(4)
		xg(4)=matr(3)
        xg=xg/sum(xg)
        close(fi)
    elseif (strval=="Qmom0D") then
        after_foaming=TRIM(ADJUSTL(qmom0D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xg(1)=0 ! no oxygen
		xg(2)=0 ! no nitrogen
        xg(3)=matr(4)
		xg(4)=matr(3)
        xg=xg/sum(xg)
        close(fi)
    elseif (strval=="Qmom3D") then
        after_foaming=TRIM(ADJUSTL(qmom3D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xg(1)=0 ! no oxygen
		xg(2)=0 ! no nitrogen
        xg(3)=matr(4)
		xg(4)=matr(3)
        xg=xg/sum(xg)
        close(fi)
    elseif (strval=="DirectInput") then
		call fson_get(json_data, "foamCondition.initialComposition.O2", xg(1))
		call fson_get(json_data, "foamCondition.initialComposition.N2", xg(2))
		call fson_get(json_data, "foamCondition.initialComposition.CO2", xg(3))
		call fson_get(&
			json_data, "foamCondition.initialComposition.Cyclopentane", xg(4))
		xg=xg/sum(xg)
    else
        write(*,*) 'unknown source for gas composition'
        stop
    endif
	call fson_get(json_data, "sourceOfProperty.foamDensity", strval)
    if (strval=="BubbleGrowth") then
        after_foaming=TRIM(ADJUSTL(bg_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        rhof=matr(1)
        close(fi)
    elseif (strval=="Qmom0D") then
        after_foaming=TRIM(ADJUSTL(qmom0D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        rhof=matr(1)
        close(fi)
    elseif (strval=="Qmom3D") then
        after_foaming=TRIM(ADJUSTL(qmom3D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        rhof=matr(1)
        close(fi)
    elseif (strval=="DirectInput") then
		call fson_get(json_data, "morphology.foamDensity", rhof)
    else
        write(*,*) 'unknown source for foam density'
        stop
    endif
	call fson_get(json_data, "sourceOfProperty.cellSize", strval)
    if (strval=="BubbleGrowth") then
        after_foaming=TRIM(ADJUSTL(bg_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        dcell=matr(2)
        close(fi)
    elseif (strval=="Qmom0D") then
        after_foaming=TRIM(ADJUSTL(qmom0D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        dcell=matr(2)
        close(fi)
    elseif (strval=="Qmom3D") then
        after_foaming=TRIM(ADJUSTL(qmom3D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        dcell=matr(2)
        close(fi)
    elseif (strval=="DirectInput") then
		call fson_get(json_data, "morphology.cellSize", dcell)
    else
        write(*,*) 'unknown source for cell size'
        stop
    endif
	call fson_get(json_data, "sourceOfProperty.strutContent", strval)
    if (strval=="StrutContent") then
        call strutContent(fstrut,rhof)
    elseif (strval=="DirectInput") then
		call fson_get(json_data, "morphology.strutContent", fstrut)
    else
        write(*,*) 'unknown source for strut content'
        stop
    endif
	call fson_get(json_data, "sourceOfProperty.wallThickness", strval)
    if (strval=="DirectInput") then
		call fson_get(json_data, "morphology.wallThickness", dwall)
    else
        write(*,*) 'unknown source for wall thickness'
        stop
    endif
	call fson_get(json_data, "physicalProperties.polymerDensity", rhop)
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.O2", strval)
	if (strval=="constant") then
		solModel(1)=0
	elseif ( strval=="modena" ) then
		solModel(1)=1
	else
		print*, "Solubility model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.N2", strval)
	if (strval=="constant") then
		solModel(2)=0
	elseif ( strval=="modena" ) then
		solModel(2)=1
	else
		print*, "Solubility model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.CO2", strval)
	if (strval=="constant") then
		solModel(3)=0
	elseif ( strval=="modena" ) then
		solModel(3)=1
	else
		print*, "Solubility model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.solubilityModel.Cyclopentane", strval)
	if (strval=="constant") then
		solModel(4)=0
	elseif ( strval=="modena" ) then
		solModel(4)=1
	else
		print*, "Solubility model must be constant or modena"
	endif
	if (solModel(1)==0) then
		call fson_get(json_data, "physicalProperties.foam.solubility.O2", Sg(1))
	endif
	if (solModel(2)==0) then
		call fson_get(json_data, "physicalProperties.foam.solubility.N2", Sg(2))
	endif
	if (solModel(3)==0) then
		call fson_get(&
			json_data, "physicalProperties.foam.solubility.CO2", Sg(3))
	endif
	if (solModel(4)==0) then
		call fson_get(&
			json_data, "physicalProperties.foam.solubility.Cyclopentane", Sg(4))
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.O2", strval)
	if (strval=="constant") then
	    diffModel(1)=0
	elseif ( strval=="modena" ) then
		diffModel(1)=1
	else
		print*, "Diffusivity model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.N2", strval)
	if (strval=="constant") then
	    diffModel(2)=0
	elseif ( strval=="modena" ) then
		diffModel(2)=1
	else
		print*, "Diffusivity model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.CO2",strval)
	if (strval=="constant") then
		diffModel(3)=0
	elseif ( strval=="modena" ) then
		diffModel(3)=1
	else
		print*, "Diffusivity model must be constant or modena"
	endif
	call fson_get(json_data, &
		"physicalProperties.foam.diffusivityModel.Cyclopentane", strval)
	if (strval=="constant") then
		diffModel(4)=0
	elseif ( strval=="modena" ) then
		diffModel(4)=1
	else
		print*, "Diffusivity model must be constant or modena"
	endif
	if (diffModel(1)==0) then
		call fson_get(json_data, &
			"physicalProperties.foam.diffusivity.O2", Dg(1))
	endif
	if (diffModel(2)==0) then
		call fson_get(json_data, &
			"physicalProperties.foam.diffusivity.N2", Dg(2))
	endif
	if (diffModel(3)==0) then
		call fson_get(json_data, &
			"physicalProperties.foam.diffusivity.CO2", Dg(3))
	endif
	if (diffModel(4)==0) then
		call fson_get(json_data, &
			"physicalProperties.foam.diffusivity.Cyclopentane", Dg(4))
	endif
	if (sheet) then
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.O2", sheetSg(1))
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.N2", sheetSg(2))
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.CO2", sheetSg(3))
		call fson_get(json_data, &
			"physicalProperties.sheet.solubility.Cyclopentane", sheetSg(4))
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.O2", sheetDg(1))
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.N2", sheetDg(2))
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.CO2", sheetDg(3))
		call fson_get(json_data, &
			"physicalProperties.sheet.diffusivity.Cyclopentane", sheetDg(4))
	endif
	call fson_destroy(json_data)
end subroutine input
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Saves results to file.
!!
!! Saves partial pressure profiles in whole foam and in gas phase only.
!! @param [in] time time
subroutine output(iprof, time, ystate, neq, keq, fi)
	use constants
	use globals
	integer, intent(in) :: iprof !< number time step
	integer, intent(in) :: neq !< number of equations
	integer, intent(in) :: fi !< file index for keq_time.out
	real(dp), intent(in) :: time !< time
	real(dp), intent(in) :: ystate(:) !< integrated variables
	real(dp), intent(in) :: keq !< equivalent conductivity
	integer :: i, j
	integer :: spp
	real(dp) :: pos
	real(dp), dimension(:), allocatable :: pp

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
	open(unit=11,file='H2perm_'//trim(name_f)//'.dat')
	open(unit=12,file='ppar_'//trim(name_f)//'.dat')

	! profiles
	do i = 1, neq/ngas
		write (11,100) time/(3600*24),pos,&
			ystate(ngas*(i-1)+1)*sol(ngas*(i-1)+1),&
			ystate(ngas*(i-1)+2)*sol(ngas*(i-1)+2),&
			ystate(ngas*(i-1)+3)*sol(ngas*(i-1)+3),&
			ystate(ngas*(i-1)+4)*sol(ngas*(i-1)+4)
	enddo
	allocate(pp(ngas))
	pp=0
	spp=0
	do i = 1,nfv
		if (i==1) then
		    pos=dz(1)
		else
			pos=pos+(dz(i-1)+dz(i))/2
		endif
		if (mor(i)==1) then
			do j=1,ngas
				pp(j)=pp(j)+ystate(ngas*(i-1)+j)*Rg*temp
			enddo
			spp=spp+1
			write (12,101) time/(3600*24),pos,ystate(ngas*(i-1)+1)*Rg*temp,&
				ystate(ngas*(i-1)+2)*Rg*temp,ystate(ngas*(i-1)+3)*Rg*temp,&
				ystate(ngas*(i-1)+4)*Rg*temp
		endif
	enddo
	pp=pp/spp
    write(fi,'(10es23.15)') tbeg/(3600*24),keq*1.0e3_dp,sum(pp),pp
    close(11)
	close(12)
100   format (f8.2,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3)
101   format (f8.2,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3,2x,ES23.8E3)
end subroutine output
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Prints header to console.
!!
!! Outputs useful information, which identifies, which foam is simulated.
subroutine print_header
	use globals
	integer :: i
	print*, 'Model type: ',TRIM(ADJUSTL(modelType))
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
	do i=1,ngas
		write(*,'(A30,EN12.3,1x,A)') &
			TRIM(ADJUSTL(gasname(i)))//' diffusivity:',Dg(i),'m2/s'
	enddo
	do i=1,ngas
		write(*,'(A30,EN12.3,1x,A)') &
			TRIM(ADJUSTL(gasname(i)))//' solubility:',Sg(i),'g/g/bar'
	enddo
	do i=1,ngas
		write(*,'(A30,EN12.3,1x,A)') &
			TRIM(ADJUSTL(gasname(i)))//' permeability:',Pg(i),'mol/m/s/Pa'
	enddo
	if (modelType == 'homogeneous') then
		do i=1,ngas
			write(*,'(A30,EN12.3,1x,A)') &
				TRIM(ADJUSTL(gasname(i)))//' effective diffusivity:',&
				Deff(i),'m2/s'
		enddo
	endif
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
