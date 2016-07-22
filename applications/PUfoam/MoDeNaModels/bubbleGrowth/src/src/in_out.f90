!> @file
!! handles input and output
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module in_out
    use foaming_globals_m
    use constants
    use globals
    use ioutils, only:newunit,str
    implicit none
    private
    public set_paths,read_inputs,save_integration_header,&
        save_integration_step,save_integration_close,load_old_results
contains
!********************************BEGINNING*************************************
!> set paths to all files
subroutine set_paths
    fileplacein='../inputs/'
    fileplaceout='./'
    inputs='unifiedInput.json'
    outputs_1d='outputs_1d.out'
    outputs_GR='outputs_GR.out'
    outputs_c='outputs_c.out'
    outputs_kin='kinetics.out'
    outputs_drain='bblgr_2_drain.out'
    inputs=TRIM(ADJUSTL(fileplacein))//TRIM(ADJUSTL(inputs))
    outputs_1d=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_1d))
    outputs_GR=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_GR))
    outputs_c=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_c))
    outputs_kin=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_kin))
    outputs_drain=TRIM(ADJUSTL(fileplaceout))//TRIM(ADJUSTL(outputs_drain))
end subroutine set_paths
!***********************************END****************************************
!> reads input values from a file
!! save them to global variables
subroutine read_inputs
    use fson
    use fson_value_m, only: fson_value_get
    character(len=1024) :: strval
    type(fson_value), pointer :: json_data
    write(*,*) 'loading input file ',TRIM(inputs)
    json_data => fson_parse(inputs)
    ngas=2
    co2_pos=2
    call fson_get(json_data, "bubbleGrowth.geometry", geometry)
    call fson_get(json_data, "bubbleGrowth.integrator", strval)
    if (strval=="dlsode") then
        integrator=1
    elseif (strval=="dlsodes") then
        integrator=2
    elseif (strval=="cvode") then
        integrator=3
    else
        stop 'unknown integrator'
    endif
    call fson_get(json_data, "bubbleGrowth.method", strval)
    if (strval=="nonstiff") then
        if (integrator==1 .or. integrator==2) then
            int_meth=10
        elseif (integrator==3) then
            int_meth=1
        endif
    elseif (strval=="stiff") then
        if (integrator==1 .or. integrator==2) then
            if (integrator==1) then
                int_meth=22
            elseif (integrator==2) then
                int_meth=222
            endif
        elseif (integrator==3) then
            int_meth=2
        endif
    else
        stop 'method can be either stiff or nonstiff'
    endif
    call fson_get(json_data, "bubbleGrowth.inertialTerm", inertial_term)
    call fson_get(json_data, "bubbleGrowth.solubilityCorrection", solcorr)
    call fson_get(json_data, "bubbleGrowth.meshCoarseningParameter", mshco)
    call fson_get(json_data, "bubbleGrowth.internalNodes", p)
    call fson_get(json_data, "bubbleGrowth.initialTime", tstart)
    if (firstrun .and. .not. shooting) call fson_get(json_data, "bubbleGrowth.finalTime", tend)
    call fson_get(json_data, "bubbleGrowth.outerTimeSteps", its)
    call fson_get(json_data, "bubbleGrowth.maxInnerTimeSteps", maxts)
    call fson_get(json_data, "bubbleGrowth.relativeTolerance", rel_tol)
    call fson_get(json_data, "bubbleGrowth.absoluteTolerance", abs_tol)
    allocate(D(ngas),cbl(ngas),xgas(ngas+1),KH(ngas),Mbl(ngas),&
        dHv(ngas),mb(ngas),mb2(ngas),mb3(ngas),avconc(ngas),pressure(ngas),&
        diff_model(ngas),sol_model(ngas),cpblg(ngas),cpbll(ngas),&
        wblpol(ngas),D0(ngas))
    call fson_get(json_data, "physicalProperties.pressure", pamb)
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.molarMass", Mbl(1))
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.molarMass", Mbl(2))
    call fson_get(json_data, "physicalProperties.polymer.heatCapacity", cppol)
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.heatCapacityInLiquidPhase", cpbll(1))
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.heatCapacityInLiquidPhase", cpbll(2))
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.heatCapacityInGaseousPhase", cpblg(1))
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.heatCapacityInGaseousPhase", cpblg(2))
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.evaporationHeat", dHv(1))
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.evaporationHeat", dHv(2))
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.density", rhobl)
    call fson_get(json_data, "initialConditions.temperature", temp0)
    if (.not. shooting) call fson_get(json_data, "initialConditions.bubbleRadius", R0)
    call fson_get(json_data, "initialConditions.numberBubbleDensity", nb0)
    call fson_get(json_data, "kinetics.kineticModel", strval)
    if (strval=='Baser') then
        kin_model=1
    elseif (strval=='BaserRx') then
        kin_model=3
    elseif (strval=='RF-1') then
        kin_model=4
    else
        stop 'unknown kinetic model'
    endif
    call fson_get(json_data, "initialConditions.concentrations.water", W0)
    call fson_get(json_data, "kinetics.gelPoint", gelpointconv)
    if (kin_model==1 .or. kin_model==3) then
        call fson_get(json_data, "kinetics.useDilution", dilution)
        call fson_get(json_data, "initialConditions.concentrations.polyol", OH0)
        call fson_get(json_data, "initialConditions.concentrations.isocyanate", NCO0)
        call fson_get(json_data, "kinetics.gellingReaction.frequentialFactor", AOH)
        call fson_get(json_data, "kinetics.gellingReaction.activationEnergy", EOH)
        call fson_get(json_data, "kinetics.blowingReaction.frequentialFactor", AW)
        call fson_get(json_data, "kinetics.blowingReaction.activationEnergy", EW)
    elseif (kin_model==4) then
        call fson_get(json_data, "initialConditions.concentrations.catalyst", catalyst)
        call fson_get(json_data, "initialConditions.concentrations.polyol1", polyol1_ini)
        call fson_get(json_data, "initialConditions.concentrations.polyol2", polyol2_ini)
        call fson_get(json_data, "initialConditions.concentrations.amine", amine_ini)
        call fson_get(json_data, "initialConditions.concentrations.isocyanate1", isocyanate1_ini)
        call fson_get(json_data, "initialConditions.concentrations.isocyanate2", isocyanate2_ini)
        call fson_get(json_data, "initialConditions.concentrations.isocyanate3", isocyanate3_ini)
        OH0=polyol1_ini+polyol2_ini
    endif
    call fson_get(json_data, "initialConditions.concentrations.blowingAgents.PBL", cbl(1))
    call fson_get(json_data, "initialConditions.concentrations.blowingAgents.CO2", cbl(2))
    if (kin_model==1 .or. kin_model==3 .or. kin_model==4) then
        call fson_get(json_data, "kinetics.gellingReaction.reactionEnthalpy", dHOH)
        call fson_get(json_data, "kinetics.blowingReaction.reactionEnthalpy", dHW)
    endif
    call fson_get(json_data, "physicalProperties.polymer.polymerDensityModel", strval)
    if (strval=="constant") then
        rhop_model=1
        call fson_get(json_data, "physicalProperties.polymer.density", rhop)
    elseif (strval=="nanotools") then
        rhop_model=2
    elseif (strval=="pcsaft") then
        rhop_model=3
    else
        stop 'unknown polymer density model'
    endif
    call fson_get(json_data, "physicalProperties.surfaceTensionModel", strval)
    if (strval=="constant") then
        itens_model=1
        call fson_get(json_data, "physicalProperties.surfaceTension", sigma)
    elseif (strval=="pcsaft") then
        itens_model=2
    else
        stop 'unknown interfacial tension model'
    endif
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.diffusivityModel", strval)
    if (strval=="constant") then
        diff_model(1)=1
        call fson_get(json_data, "physicalProperties.blowingAgents.PBL.diffusivity", D(1))
    elseif (strval=="nanotools") then
        diff_model(1)=2
    else
        stop 'unknown diffusivity model'
    endif
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.diffusivityModel", strval)
    if (strval=="constant") then
        diff_model(2)=1
        call fson_get(json_data, "physicalProperties.blowingAgents.CO2.diffusivity", D(2))
    elseif (strval=="nanotools") then
        diff_model(2)=2
    else
        stop 'unknown diffusivity model'
    endif
    call fson_get(json_data, "physicalProperties.blowingAgents.PBL.solubilityModel", strval)
    if (strval=="constant") then
        sol_model(1)=1
        call fson_get(json_data, "physicalProperties.blowingAgents.PBL.solubility", KH(1))
    elseif (strval=="pcsaft") then
        sol_model(1)=2
    elseif (strval=="Gupta") then
        sol_model(1)=3
    elseif (strval=="Winkler") then
        sol_model(1)=4
    elseif (strval=="Baser") then
        sol_model(1)=6
    elseif (strval=="hardcodedWinkler") then
        sol_model(1)=7
    else
        stop 'unknown solubility model'
    endif
    call fson_get(json_data, "physicalProperties.blowingAgents.CO2.solubilityModel", strval)
    if (strval=="constant") then
        sol_model(2)=1
    elseif (strval=="pcsaft") then
        sol_model(2)=2
    elseif (strval=="hardcodedconstant") then
        sol_model(2)=8
        call fson_get(json_data, "physicalProperties.blowingAgents.CO2.solubility", KH(2))
    else
        stop 'unknown solubility model'
    endif
    call fson_get(json_data, "physicalProperties.polymer.viscosityModel", strval)
    if (strval=="constant") then
        visc_model=1
        call fson_get(json_data, "physicalProperties.polymer.viscosity", eta)
    elseif (strval=="hardcodedCastroMacosko") then
        visc_model=2
    elseif (strval=="CastroMacosko") then
        visc_model=3
    else
        stop 'unknown viscosity model'
    endif
    call fson_destroy(json_data)
    write(*,*) 'done: inputs loaded'
    write(*,*)
end subroutine read_inputs
!***********************************END****************************************


!********************************BEGINNING*************************************
!> opens output files and writes a header
subroutine save_integration_header
    open (unit=newunit(fi1), file = outputs_1d)
    write(fi1,'(1000A24)') '#time', 'radius','pressure1', 'pressure2',&
        'conversion_of_polyol',&
        'conversion_of_water', 'eq_concentration', 'first_concentration', &
        'viscosity', 'moles_in_polymer', 'moles_in_bubble', 'total_moles', &
        'shell_thickness', 'temperature', 'foam_density', 'weight_fraction1', &
        'weight_fraction2','porosity'
    open (unit=newunit(fi2), file = outputs_GR)
    write(fi2,'(1000A23)') '#GrowthRate1', 'GrowthRate2', 'temperature', &
        'bubbleRadius','c1','c2','p1','p2'
    if (kin_model==4) then
        open (unit=newunit(fi3), file = outputs_kin)
        write(fi3,'(1000A23)') "time","Catalyst_1","CE_A0","CE_A1","CE_B",&
            "CE_B2","CE_I0","CE_I1","CE_I2","CE_PBA","CE_Breac","CE_Areac0",&
            "CE_Areac1","CE_Ireac0","CE_Ireac1","CE_Ireac2","Bulk","R_1",&
            "R_1_mass","R_1_temp","R_1_vol"
    endif
    open (unit=newunit(fi4), file = outputs_c)
    open (unit=newunit(fi5), file = outputs_drain)
    write(fi5,'(1000A24)') '#time','radius','porosity','viscosity'
end subroutine save_integration_header
!***********************************END****************************************


!********************************BEGINNING*************************************
!> writes an integration step to output file
subroutine save_integration_step(iout)
    integer :: i,iout
    write(fi1,"(1000es24.15e3)") time,radius,pressure,Y(xOHeq),Y(xWeq),&
        eqconc,Y(fceq),eta,mb(1),mb2(1),mb3(1),st,temp,rhofoam,&
        wblpol,porosity
    write(fi2,"(1000es23.15)") grrate, temp, radius, avconc, pressure
    if (kin_model==4) then
        write(fi3,"(1000es23.15)") time,Y(kineq(1)),Y(kineq(2)),Y(kineq(3)),&
            Y(kineq(4)),Y(kineq(5)),Y(kineq(6)),Y(kineq(7)),Y(kineq(8)),&
            Y(kineq(9)),Y(kineq(10)),Y(kineq(11)),Y(kineq(12)),Y(kineq(13)),&
            Y(kineq(14)),Y(kineq(15)),Y(kineq(16)),Y(kineq(17)),Y(kineq(18)),&
            Y(kineq(19)),Y(kineq(20))
    endif
    write(fi4,"(1000es23.15)") (Y(fceq+i+1),i=0,ngas*p,ngas)
    write(fi5,"(1000es24.15e3)") time,radius,porosity,eta
    ! save arrays, which are preserved for future use
    etat(iout,1)=time
    etat(iout,2)=eta
    port(iout,1)=time
    port(iout,2)=porosity
    init_bub_rad(iout,1)=time
    init_bub_rad(iout,2)=radius
end subroutine save_integration_step
!***********************************END****************************************


!********************************BEGINNING*************************************
!> closes output files
subroutine save_integration_close(iout)
    integer :: iout
    real(dp), dimension(:,:), allocatable :: matr
    close(fi1)
    close(fi2)
    if (kin_model==4) then
        close(fi3)
    endif
    close(fi4)
    close(fi5)
    ! reallocate matrices for eta_rm and bub_vf functions
    ! interpolation doesn't work otherwise
    if (iout /= its) then
        allocate(matr(0:its,2))
        matr=etat
        deallocate(etat)
        allocate(etat(0:iout,2))
        etat=matr(0:iout,:)
        matr=port
        deallocate(port)
        allocate(port(0:iout,2))
        port=matr(0:iout,:)
        matr=init_bub_rad
        deallocate(init_bub_rad)
        allocate(init_bub_rad(0:iout,2))
        init_bub_rad=matr(0:iout,:)
    endif
end subroutine save_integration_close
!***********************************END****************************************


!********************************BEGINNING*************************************
!> loads old results
subroutine load_old_results
    integer :: i,j,ios,fi
    real(dp), dimension(:,:), allocatable :: matrix
    j=0
    open(newunit(fi),file=outputs_1d)
        do  !find number of points
            read(fi,*,iostat=ios)
            if (ios/=0) exit
            j=j+1
        enddo
        allocate(matrix(j,18))
        rewind(fi)
        read(fi,*)
        do i=2,j
            read(fi,*) matrix(i,:)
        enddo
    close(fi)
    allocate(bub_rad(j-1,2))
    bub_rad(:,1)=matrix(2:j,1)
    bub_rad(:,2)=matrix(2:j,2)
    deallocate(matrix)
end subroutine load_old_results
!***********************************END****************************************
end module in_out
