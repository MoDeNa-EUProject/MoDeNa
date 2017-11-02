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
	public input,output,print_header,write_header
contains
!********************************BEGINNING*************************************
!> Reads input file.
!!
!! Inputs are saved to global variables.
subroutine input()
	use constants
	use globals
	use physicalProperties
	use ioutils, only: newunit
	use fson
    use fson_value_m, only: fson_value_get
	type(fson_value), pointer :: json_data
	character(len=1024) :: strval
    character(len=99) :: after_foaming,after_foaming0='after_foaming.txt'
    character(len=99) :: bg_res='../../foamExpansion/results/bubbleGrowth/'
    character(len=99) :: qmom0D_res='../../foamExpansion/results/CFD0D/'
    character(len=99) :: qmom3D_res='../../foamExpansion/results/CFD3D/'
	real(dp) :: matr(7)
	integer :: fi, i
	call get_names(gasname)
	ngas=size(gasname)
    allocate(solModel(ngas),diffModel(ngas))
    allocate(Sg(ngas),Dg(ngas),Pg(ngas),Deff(ngas),Seff(ngas))
    allocate(sheetSg(ngas),sheetDg(ngas))
    allocate(pBg(ngas),xg(ngas),kfoamXg(ngas),kgasXg(ngas))
    allocate(sgModena(ngas),sgInputs(ngas),sgOutputs(ngas))
    allocate(sgTemppos(ngas),sgxl1pos(ngas),sgxl2pos(ngas))
    allocate(dgModena(ngas),dgInputs(ngas),dgOutputs(ngas),dgTemppos(ngas))
    allocate(kgModena(ngas),kgInputs(ngas),kgOutputs(ngas),kgTemppos(ngas))
	allocate(Mg(ngas))
	! Read input parameters
	json_data => fson_parse("../inputs/foamAging.json")
	call fson_get(json_data, "modelType", modelType)
	if (modelType /= "heterogeneous" .and. modelType /= "homogeneous") then
		write(*,*) 'modelType must be heterogeneous or homogeneous'
		stop
	endif
	call fson_get(json_data, "foamCondition.inProtectiveSheet", sheet)
	call fson_get(json_data, "numerics.timeStart", tbeg)
	call fson_get(json_data, "numerics.timeEnd", tend)
	call fson_get(json_data, "numerics.progressTime", progressTime)
	if (progressTime == "linear") then
		call fson_get(json_data, "numerics.numberOfOutputs", nroutputs)
    elseif (progressTime == "logarithmic") then
		call fson_get(json_data, "numerics.outputsPerOrder", outputsPerOrder)
		call fson_get(json_data, "numerics.numberOfOrders", numberOfOrders)
        nroutputs = numberOfOrders*outputsPerOrder + 1
	else
		write(*,*) 'modelType must be linear or logarithmic'
		stop
    endif
	if (modelType == "heterogeneous") then
		call fson_get(json_data, "numerics.wallDiscretization", divwall)
		call fson_get(json_data, "numerics.cellDiscretization", divcell)
	elseif (modelType == "homogeneous") then
		call fson_get(json_data, "numerics.foamDiscretization", divfoam)
	endif
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
	do i = 1, ngas
		call fson_get(json_data, "foamCondition.boundaryPressure."//&
			trim(adjustl(gasname(i))), pBg(i))
	enddo
	call fson_get(json_data, "sourceOfProperty.gasComposition", strval)
    if (strval=="BubbleGrowth") then
        after_foaming=TRIM(ADJUSTL(bg_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xg = 0
        ! this works only for cyclopentane and CO2
        xg(1) = matr(4)
        xg(2) = matr(3)
        close(fi)
    elseif (strval=="Qmom0D") then
        after_foaming=TRIM(ADJUSTL(qmom0D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xg = 0
        ! this works only for cyclopentane and CO2
        xg(1) = matr(4)
        xg(2) = matr(3)
        close(fi)
    elseif (strval=="Qmom3D") then
        after_foaming=TRIM(ADJUSTL(qmom3D_res))//TRIM(ADJUSTL(after_foaming0))
        open(unit=newunit(fi),file=after_foaming)
        read(fi,*)
        read(fi,*) matr(1:4)
        xg = 0
        ! this works only for cyclopentane and CO2
        xg(1) = matr(4)
        xg(2) = matr(3)
        close(fi)
    elseif (strval=="DirectInput") then
		do i=1,ngas
            call fson_get(&
                json_data,&
                "foamCondition.initialComposition."//TRIM(ADJUSTL(gasname(i))),&
                xg(i))
        enddo
    else
        write(*,*) 'unknown source for gas composition'
        stop
    endif
	xg = xg / sum(xg)
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
	do i = 1, ngas
		call fson_get(json_data, &
			"physicalProperties.molarMass."//&
			trim(adjustl(gasname(i))), Mg(i))
		call fson_get(json_data, &
			"physicalProperties.foam.solubilityModel."//&
			trim(adjustl(gasname(i))), strval)
		if (strval=="constant") then
			solModel(i)=0
			call fson_get(json_data, &
				"physicalProperties.foam.solubility."//&
				trim(adjustl(gasname(i))), Sg(i))
		elseif ( strval=="modena" ) then
			solModel(i)=1
		else
			print*, "Solubility model must be constant or modena"
		endif
		call fson_get(json_data, &
			"physicalProperties.foam.diffusivityModel."//&
			trim(adjustl(gasname(i))), strval)
		if (strval=="constant") then
			diffModel(i)=0
			call fson_get(json_data, &
				"physicalProperties.foam.diffusivity."//&
			trim(adjustl(gasname(i))), Dg(i))
		elseif ( strval=="modena" ) then
			diffModel(i)=1
		elseif (strval == "foam") then
			diffModel(i) = 2
			call fson_get(json_data, &
				"physicalProperties.foam.diffusivity."//&
				trim(adjustl(gasname(i))), Deff(i))
		else
			print*, "Diffusivity model must be constant, foam or modena"
		endif
		if (sheet) then
			call fson_get(json_data, &
				"physicalProperties.sheet.solubility."//&
				trim(adjustl(gasname(i))), sheetSg(i))
			call fson_get(json_data, &
				"physicalProperties.sheet.diffusivity."//&
				trim(adjustl(gasname(i))), sheetDg(i))
		endif
	enddo
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
	use ioutils, only: newunit
	integer, intent(in) :: iprof !< number time step
	integer, intent(in) :: neq !< number of equations
	integer, intent(in) :: fi !< file index for degas_scalar.csv
	real(dp), intent(in) :: time !< time
	real(dp), intent(in) :: ystate(:) !< integrated variables
	real(dp), intent(in) :: keq !< equivalent conductivity
	integer :: i, j, fi2, fi3
	integer :: spp
	real(dp) :: pos
	real(dp), dimension(:), allocatable :: pp,vals
	character(len=255) :: fmt ! format for writing

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

	allocate(vals(ngas))
	if (modelType == "heterogeneous") then
		open(newunit(fi2),file='conc_'//trim(name_f)//'.csv')
		write(fmt,*) ngas + 1
		fmt = '(1x, A23, '//trim(adjustl(fmt))//'(",", A23))'
		write (fi2,fmt) 'time', 'position', gasname
		write(fmt,*) ngas + 1
		fmt = '(1x, ES23.15, '//trim(adjustl(fmt))//'(",", ES23.15))'
		do i = 1, neq/ngas
			do j = 1, ngas
				vals(j) = ystate(ngas*(i-1)+j)*sol(ngas*(i-1)+j)
			enddo
			write (fi2,fmt) time/(3600*24),pos,vals
		enddo
    	close(fi2)
	endif

	open(newunit(fi3),file='pres_'//trim(name_f)//'.csv')
	write(fmt,*) ngas + 1
	fmt = '(1x, A23, '//trim(adjustl(fmt))//'(",", A23))'
	write (fi3,fmt) 'time', 'position', gasname
	write(fmt,*) ngas + 1
	fmt = '(1x, ES23.15, '//trim(adjustl(fmt))//'(",", ES23.15))'
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
				if (modelType == "heterogeneous") then
					pp(j)=pp(j)+ystate(ngas*(i-1)+j)*Rg*temp
				elseif ( modelType == "homogeneous" ) then
					pp(j)=pp(j)+ystate(ngas*(i-1)+j)/Seff(j)
				endif
			enddo
			spp=spp+1
			if (modelType == "heterogeneous") then
				do j = 1, ngas
					vals(j) = ystate(ngas*(i-1)+j)*Rg*temp
				enddo
				write (fi3,fmt) time/(3600*24),pos,vals
			elseif ( modelType == "homogeneous" ) then
				do j = 1, ngas
					vals(j) = ystate(ngas*(i-1)+j)/Seff(j)
				enddo
				write (fi3,fmt) &
					time/(3600*24),pos,vals
			endif
		endif
	enddo
	close(fi3)
	pp=pp/spp
	write(fmt,*) ngas + 2
	fmt = '(1x, ES23.15, '//trim(adjustl(fmt))//'(",", ES23.15))'
    write(fi,fmt) time/(3600*24),keq*1.0e3_dp,sum(pp),pp
	deallocate(vals)
end subroutine output
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Writes header of degas_scalar.csv.
!!
!! Outputs evolution of conductivity and average gas composition.
subroutine write_header(fi)
	use globals
	use ioutils, only: newunit
	integer, intent(out) :: fi
	character(len=255), dimension(:), allocatable :: pres_names(:)
	integer :: i
	character(len=255) :: fmt
	allocate(pres_names(ngas))
	do i = 1, ngas
		pres_names(i) = 'p_'//trim(adjustl(gasname(i)))
	enddo
	write(fmt,*) ngas + 2
	fmt = '(1x, A23, '//trim(adjustl(fmt))//'(",", A23))'
	open (newunit(fi),file='degas_scalar.csv')
    write(fi, fmt) 'time', 'eq_conductivity', 'total_pressure', pres_names
	deallocate(pres_names)
end subroutine write_header
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
				TRIM(ADJUSTL(gasname(i)))//' effective solubility:',&
				Seff(i),'mol/m3/Pa'
		enddo
		do i=1,ngas
			write(*,'(A30,EN12.3,1x,A)') &
				TRIM(ADJUSTL(gasname(i)))//' effective diffusivity:',&
				Deff(i),'m2/s'
		enddo
	endif
    print*, 'Numerics:'
	if (modelType == "heterogeneous") then
	    write(*,'(A30,I12)') 'finite volumes in wall:',divwall
	    write(*,'(A30,I12)') 'finite volumes in cell:',divcell
	elseif ( modelType == "homogeneous" ) then
		write(*,'(A30,I12)') 'finite volumes in foam:',divfoam
	endif
    if (sheet) then
        write(*,'(A30,I12)') 'finite volumes in sheet:',divsheet
    endif
    write(*,'(A30,I12)') 'finite volumes in total:',nfv
    write(*,'(A30,EN12.3,1x,A)') 'initial time:',tbeg,'s'
    write(*,'(A30,EN12.3,1x,A)') 'end time:',tend,'s'
end subroutine print_header
!***********************************END****************************************
end module inout
