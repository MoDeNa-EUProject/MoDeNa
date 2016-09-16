!> @file
!! subroutines for calculation of physical properties of polymer and blowing
!! agents using Modena calls
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module physicalProperties
    use constants
    use globals, only: solModel,diffModel
    use fmodena
    implicit none
    !modena variables
    integer(c_int) :: ret
    type(c_ptr) :: rhopModena = c_null_ptr
    type(c_ptr) :: rhopInputs = c_null_ptr
    type(c_ptr) :: rhopOutputs = c_null_ptr
    integer(c_size_t) :: rhopTemppos
    type(c_ptr) :: kfoamModena = c_null_ptr
    type(c_ptr) :: kfoamInputs = c_null_ptr
    type(c_ptr) :: kfoamOutputs = c_null_ptr
    integer(c_size_t) :: kfoamEpspos
    integer(c_size_t) :: kfoamDcellpos
    integer(c_size_t) :: kfoamFstrutpos
    ! integer(c_size_t) :: kfoamKgaspos
    integer(c_size_t) :: kfoamTemppos
    integer(c_size_t) :: kfoamXCO2pos
    integer(c_size_t) :: kfoamXCyPpos
    integer(c_size_t) :: kfoamXO2pos
    integer(c_size_t) :: kfoamXN2pos
    type(c_ptr) :: kgasModena = c_null_ptr
    type(c_ptr) :: kgasInputs = c_null_ptr
    type(c_ptr) :: kgasOutputs = c_null_ptr
    integer(c_size_t) :: kgasTemppos
    integer(c_size_t) :: kgasXCO2pos
    integer(c_size_t) :: kgasXCyPpos
    integer(c_size_t) :: kgasXO2pos
    integer(c_size_t) :: kgasXN2pos
    type(c_ptr) :: kcdModena = c_null_ptr
    type(c_ptr) :: kcdInputs = c_null_ptr
    type(c_ptr) :: kcdOutputs = c_null_ptr
    integer(c_size_t) :: kcdTemppos
    type(c_ptr) :: kairModena = c_null_ptr
    type(c_ptr) :: kairInputs = c_null_ptr
    type(c_ptr) :: kairOutputs = c_null_ptr
    integer(c_size_t) :: kairTemppos
    type(c_ptr) :: kcypModena = c_null_ptr
    type(c_ptr) :: kcypInputs = c_null_ptr
    type(c_ptr) :: kcypOutputs = c_null_ptr
    integer(c_size_t) :: kcypTemppos
    type(c_ptr) :: scdModena = c_null_ptr
    type(c_ptr) :: scdInputs = c_null_ptr
    type(c_ptr) :: scdOutputs = c_null_ptr
    integer(c_size_t) :: scdTemppos
    integer(c_size_t) :: scdxl1pos
    integer(c_size_t) :: scdxl2pos
    type(c_ptr) :: sairModena = c_null_ptr
    type(c_ptr) :: sairInputs = c_null_ptr
    type(c_ptr) :: sairOutputs = c_null_ptr
    integer(c_size_t) :: sairTemppos
    integer(c_size_t) :: sairxl1pos
    integer(c_size_t) :: sairxl2pos
    type(c_ptr) :: scypModena = c_null_ptr
    type(c_ptr) :: scypInputs = c_null_ptr
    type(c_ptr) :: scypOutputs = c_null_ptr
    integer(c_size_t) :: scypTemppos
    integer(c_size_t) :: scypxl1pos
    integer(c_size_t) :: scypxl2pos
    type(c_ptr) :: dcdModena = c_null_ptr
    type(c_ptr) :: dcdInputs = c_null_ptr
    type(c_ptr) :: dcdOutputs = c_null_ptr
    integer(c_size_t) :: dcdTemppos
    type(c_ptr) :: dcypModena = c_null_ptr
    type(c_ptr) :: dcypInputs = c_null_ptr
    type(c_ptr) :: dcypOutputs = c_null_ptr
    integer(c_size_t) :: dcypTemppos
    type(c_ptr) :: do2Modena = c_null_ptr
    type(c_ptr) :: do2Inputs = c_null_ptr
    type(c_ptr) :: do2Outputs = c_null_ptr
    integer(c_size_t) :: do2Temppos
    type(c_ptr) :: dn2Modena = c_null_ptr
    type(c_ptr) :: dn2Inputs = c_null_ptr
    type(c_ptr) :: dn2Outputs = c_null_ptr
    integer(c_size_t) :: dn2Temppos
contains
!********************************BEGINNING*************************************
!> creates Modena models
subroutine createModels
    kfoamModena = modena_model_new (c_char_"foamConductivity"//c_null_char);
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    kfoamInputs = modena_inputs_new (kfoamModena);
    kfoamOutputs = modena_outputs_new (kfoamModena);
    kfoamEpspos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"eps"//c_null_char);
    kfoamDcellpos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"dcell"//c_null_char);
    kfoamFstrutpos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"fstrut"//c_null_char);
    ! kfoamKgaspos = modena_model_inputs_argPos(&
    !     kfoamModena, c_char_"kgas"//c_null_char);
    kfoamTemppos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"T"//c_null_char);
    kfoamXCO2pos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"x[CO2]"//c_null_char);
    kfoamXCyPpos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"x[CyP]"//c_null_char);
    kfoamXO2pos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"x[O2]"//c_null_char);
    kfoamXN2pos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"x[N2]"//c_null_char);
    call modena_model_argPos_check(kfoamModena)
    kgasModena = modena_model_new (&
        c_char_"gasMixtureConductivity"//c_null_char)
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    kgasInputs = modena_inputs_new (kgasModena)
    kgasOutputs = modena_outputs_new (kgasModena)
    kgasTemppos = modena_model_inputs_argPos(&
        kgasModena, c_char_"T"//c_null_char)
    kgasXCO2pos = modena_model_inputs_argPos(&
        kgasModena, c_char_"x[A=CO2]"//c_null_char)
    kgasXCyPpos = modena_model_inputs_argPos(&
        kgasModena, c_char_"x[A=CyP]"//c_null_char)
    kgasXO2pos = modena_model_inputs_argPos(&
        kgasModena, c_char_"x[A=O2]"//c_null_char)
    kgasXN2pos = modena_model_inputs_argPos(&
        kgasModena, c_char_"x[A=N2]"//c_null_char)
    call modena_model_argPos_check(kgasModena)
    kcdModena = modena_model_new (&
        c_char_"gas_thermal_conductivity[A=CO2]"//c_null_char);
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    kcdInputs = modena_inputs_new (kcdModena);
    kcdOutputs = modena_outputs_new (kcdModena);
    kcdTemppos = modena_model_inputs_argPos(&
        kcdModena, c_char_"T"//c_null_char);
    call modena_model_argPos_check(kcdModena)
    kairModena = modena_model_new (&
        c_char_"gas_thermal_conductivity[A=Air]"//c_null_char);
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    kairInputs = modena_inputs_new (kairModena);
    kairOutputs = modena_outputs_new (kairModena);
    kairTemppos = modena_model_inputs_argPos(&
        kairModena, c_char_"T"//c_null_char);
    call modena_model_argPos_check(kairModena)
    kcypModena = modena_model_new (&
        c_char_"gas_thermal_conductivity[A=CyP]"//c_null_char);
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    kcypInputs = modena_inputs_new (kcypModena);
    kcypOutputs = modena_outputs_new (kcypModena);
    kcypTemppos = modena_model_inputs_argPos(&
        kcypModena, c_char_"T"//c_null_char);
    call modena_model_argPos_check(kcypModena)
    if (solModel(1)==1) then
        sairModena = modena_model_new (&
            c_char_"Solubility[A=Air,B=2]"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        sairInputs = modena_inputs_new (sairModena);
        sairOutputs = modena_outputs_new (sairModena);
        sairTemppos = modena_model_inputs_argPos(&
            sairModena, c_char_"T"//c_null_char);
        sairxl1pos = modena_model_inputs_argPos(&
            sairModena, c_char_"xl1"//c_null_char);
        sairxl2pos = modena_model_inputs_argPos(&
            sairModena, c_char_"xl2"//c_null_char);
        call modena_model_argPos_check(sairModena)
    endif
    if (solModel(2)==1) then
        scdModena = modena_model_new (&
            c_char_"Solubility[A=CO2,B=2]"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        scdInputs = modena_inputs_new (scdModena);
        scdOutputs = modena_outputs_new (scdModena);
        scdTemppos = modena_model_inputs_argPos(&
            scdModena, c_char_"T"//c_null_char);
        scdxl1pos = modena_model_inputs_argPos(&
            scdModena, c_char_"xl1"//c_null_char);
        scdxl2pos = modena_model_inputs_argPos(&
            scdModena, c_char_"xl2"//c_null_char);
        call modena_model_argPos_check(scdModena)
    endif
    if (solModel(3)==1) then
        scypModena = modena_model_new (&
            c_char_"Solubility[A=CyP,B=2]"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        scypInputs = modena_inputs_new (scypModena);
        scypOutputs = modena_outputs_new (scypModena);
        scypTemppos = modena_model_inputs_argPos(&
            scypModena, c_char_"T"//c_null_char);
        scypxl1pos = modena_model_inputs_argPos(&
            scypModena, c_char_"xl1"//c_null_char);
        scypxl2pos = modena_model_inputs_argPos(&
            scypModena, c_char_"xl2"//c_null_char);
        call modena_model_argPos_check(scypModena)
    endif
    if (diffModel(1)==1) then
        do2Modena = modena_model_new (&
            c_char_"diffusivityPol[A=O2]"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        do2Inputs = modena_inputs_new (do2Modena);
        do2Outputs = modena_outputs_new (do2Modena);
        do2Temppos = modena_model_inputs_argPos(&
            do2Modena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(do2Modena)
        dn2Modena = modena_model_new (&
            c_char_"diffusivityPol[A=N2]"//c_null_char);
        dn2Inputs = modena_inputs_new (dn2Modena);
        dn2Outputs = modena_outputs_new (dn2Modena);
        dn2Temppos = modena_model_inputs_argPos(&
            dn2Modena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(dn2Modena)
    endif
    if (diffModel(2)==1) then
        dcdModena = modena_model_new (&
            c_char_"diffusivityPol[A=CO2]"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        dcdInputs = modena_inputs_new (dcdModena);
        dcdOutputs = modena_outputs_new (dcdModena);
        dcdTemppos = modena_model_inputs_argPos(&
            dcdModena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(dcdModena)
    endif
    if (diffModel(3)==1) then
        dcypModena = modena_model_new (&
            c_char_"diffusivityPol[A=CyP]"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        dcypInputs = modena_inputs_new (dcypModena);
        dcypOutputs = modena_outputs_new (dcypModena);
        dcypTemppos = modena_model_inputs_argPos(&
            dcypModena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(dcypModena)
    endif
end subroutine createModels
!***********************************END****************************************


!********************************BEGINNING*************************************
!> destroys Modena models
subroutine destroyModels
    call modena_inputs_destroy (kfoamInputs);
    call modena_outputs_destroy (kfoamOutputs);
    call modena_model_destroy (kfoamModena);
    call modena_inputs_destroy (kgasInputs);
    call modena_outputs_destroy (kgasOutputs);
    call modena_model_destroy (kgasModena);
    call modena_inputs_destroy (kcdInputs);
    call modena_outputs_destroy (kcdOutputs);
    call modena_model_destroy (kcdModena);
    call modena_inputs_destroy (kairInputs);
    call modena_outputs_destroy (kairOutputs);
    call modena_model_destroy (kairModena);
    call modena_inputs_destroy (kcypInputs);
    call modena_outputs_destroy (kcypOutputs);
    call modena_model_destroy (kcypModena);
    if (solModel(1)==1) then
        call modena_inputs_destroy (sairInputs);
        call modena_outputs_destroy (sairOutputs);
        call modena_model_destroy (sairModena);
    endif
    if (solModel(2)==1) then
        call modena_inputs_destroy (scdInputs);
        call modena_outputs_destroy (scdOutputs);
        call modena_model_destroy (scdModena);
    endif
    if (solModel(3)==1) then
        call modena_inputs_destroy (scypInputs);
        call modena_outputs_destroy (scypOutputs);
        call modena_model_destroy (scypModena);
    endif
    if (diffModel(1)==1) then
        call modena_inputs_destroy (do2Inputs);
        call modena_outputs_destroy (do2Outputs);
        call modena_model_destroy (do2Modena);
        call modena_inputs_destroy (dn2Inputs);
        call modena_outputs_destroy (dn2Outputs);
        call modena_model_destroy (dn2Modena);
    endif
    if (diffModel(2)==1) then
        call modena_inputs_destroy (dcdInputs);
        call modena_outputs_destroy (dcdOutputs);
        call modena_model_destroy (dcdModena);
    endif
    if (diffModel(3)==1) then
        call modena_inputs_destroy (dcypInputs);
        call modena_outputs_destroy (dcypOutputs);
        call modena_model_destroy (dcypModena);
    endif
end subroutine destroyModels
!***********************************END****************************************


!********************************BEGINNING*************************************
!> calculation of density of polymer
real(dp) function polymerDensity(temp)
    real(dp), intent(in) :: temp
    call modena_inputs_set(rhopInputs, rhopTemppos, temp)
    ret = modena_model_call (rhopModena, rhopInputs, rhopOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    polymerDensity=modena_outputs_get(rhopOutputs, 0_c_size_t)
end function polymerDensity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> thermal conductivity of mixture of blowing agents
real(dp) function gasConductivity(temp,xco2,xair,xcyp)
    real(dp), intent(in) :: temp,xco2,xair,xcyp
    call modena_inputs_set(kgasInputs, kgasTemppos, temp)
    call modena_inputs_set(kgasInputs, kgasXCO2pos, xCO2)
    call modena_inputs_set(kgasInputs, kgasXCyPpos, xCyP)
    call modena_inputs_set(kgasInputs, kgasXO2pos, xAir*0.21_dp)
    call modena_inputs_set(kgasInputs, kgasXN2pos, xAir*0.79_dp)
    ret = modena_model_call (kgasModena, kgasInputs, kgasOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    gasConductivity=modena_outputs_get(kgasOutputs, 0_c_size_t)
end function gasConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> thermal conductivity of carbon dioxide
real(dp) function cdConductivity(temp)
    real(dp), intent(in) :: temp
    call modena_inputs_set(kcdInputs, kcdTemppos, temp)
    ret = modena_model_call (kcdModena, kcdInputs, kcdOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    cdConductivity=modena_outputs_get(kcdOutputs, 0_c_size_t)
end function cdConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> thermal conductivity of air
real(dp) function airConductivity(temp)
    real(dp), intent(in) :: temp
    call modena_inputs_set(kairInputs, kairTemppos, temp)
    ret = modena_model_call (kairModena, kairInputs, kairOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    airConductivity=modena_outputs_get(kairOutputs, 0_c_size_t)
end function airConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> thermal conductivity of cyclo pentane
real(dp) function cypConductivity(temp)
    real(dp), intent(in) :: temp
    call modena_inputs_set(kcypInputs, kcypTemppos, temp)
    ret = modena_model_call (kcypModena, kcypInputs, kcypOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    cypConductivity=modena_outputs_get(kcypOutputs, 0_c_size_t)
end function cypConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> thermal conductivity of nitrogen
real(dp) function nitrConductivity(temp)
    real(dp), intent(in) :: temp
    nitrConductivity=(0.3918e-3_dp)+(0.9814e-4_dp)*temp+&
        (-5.0660e-8_dp)*temp**2+(1.503479e-11_dp)*temp**3
end function nitrConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> thermal conductivity of oxygen
real(dp) function oxyConductivity(temp)
    real(dp), intent(in) :: temp
    oxyConductivity=(-0.3272e-3_dp)+(0.9965e-4_dp)*temp+&
        (-3.7426e-8_dp)*temp**2+(0.973012e-11_dp)*temp**3
end function oxyConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> solubility of carbon dioxide
real(dp) function cdSolubility(temp)
    real(dp), intent(in) :: temp
    real(dp) :: xl1,xl2
    xl1=1.0e-3_dp
    xl2=1-xl1
    call modena_inputs_set(scdInputs, scdTemppos, temp)
    call modena_inputs_set(scdInputs, scdxl1pos, xl1)
    call modena_inputs_set(scdInputs, scdxl2pos, xl2)
    ret = modena_model_call (scdModena, scdInputs, scdOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    cdSolubility=modena_outputs_get(scdOutputs, 0_c_size_t)
end function cdSolubility
!***********************************END****************************************


!********************************BEGINNING*************************************
!> solubility of air
real(dp) function airSolubility(temp)
    real(dp), intent(in) :: temp
    real(dp) :: xl1,xl2
    xl1=1.0e-3_dp
    xl2=1-xl1
    call modena_inputs_set(sairInputs, sairTemppos, temp)
    call modena_inputs_set(sairInputs, sairxl1pos, xl1)
    call modena_inputs_set(sairInputs, sairxl2pos, xl2)
    ret = modena_model_call (sairModena, sairInputs, sairOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    airSolubility=modena_outputs_get(sairOutputs, 0_c_size_t)
end function airSolubility
!***********************************END****************************************


!********************************BEGINNING*************************************
!> solubility of cyclo pentane
real(dp) function cypSolubility(temp)
    real(dp), intent(in) :: temp
    real(dp) :: xl1,xl2
    xl1=1.0e-3_dp
    xl2=1-xl1
    call modena_inputs_set(scypInputs, scypTemppos, temp)
    call modena_inputs_set(scypInputs, scypxl1pos, xl1)
    call modena_inputs_set(scypInputs, scypxl2pos, xl2)
    ret = modena_model_call (scypModena, scypInputs, scypOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    cypSolubility=modena_outputs_get(scypOutputs, 0_c_size_t)
end function cypSolubility
!***********************************END****************************************


!********************************BEGINNING*************************************
!> diffusivity of carbon dioxide
real(dp) function cdDiffusivity(temp)
    real(dp), intent(in) :: temp
    call modena_inputs_set(dcdInputs, dcdTemppos, temp)
    ret = modena_model_call (dcdModena, dcdInputs, dcdOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    cdDiffusivity=modena_outputs_get(dcdOutputs, 0_c_size_t)
end function cdDiffusivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> diffusivity of cyclo pentane
real(dp) function cypDiffusivity(temp)
    real(dp), intent(in) :: temp
    call modena_inputs_set(dcypInputs, dcypTemppos, temp)
    ret = modena_model_call (dcypModena, dcypInputs, dcypOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    cypDiffusivity=modena_outputs_get(dcypOutputs, 0_c_size_t)
end function cypDiffusivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> diffusivity of air
real(dp) function airDiffusivity(temp)
    real(dp), intent(in) :: temp
    call modena_inputs_set(do2Inputs, do2Temppos, temp)
    ret = modena_model_call (do2Modena, do2Inputs, do2Outputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    airDiffusivity=0.21_dp*modena_outputs_get(do2Outputs, 0_c_size_t)
    call modena_inputs_set(dn2Inputs, dn2Temppos, temp)
    ret = modena_model_call (dn2Modena, dn2Inputs, dn2Outputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    airDiffusivity=&
        airDiffusivity+0.79_dp*modena_outputs_get(dn2Outputs, 0_c_size_t)
end function airDiffusivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> diffusivity of gases in gas phase
!! accoriding to Bird 1975, p.505, eq. 16.3-1
real(dp) function gasDiffusivity(temp)
    use globals, only: pressure
    real(dp), intent(in) :: temp
	real(dp) :: pcA, pcB, pcApcB, TcA, TcB, TcATcB
	real(dp) :: MA, MB, Mterm ,a, b, aToverTcsb
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
    gasDiffusivity = (aToverTcsb*pcApcB*TcATcB*Mterm)*1.0e5_dp/pressure
    gasDiffusivity = gasDiffusivity * 1.0e-4_dp ! m2/s
end function gasDiffusivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> heat capacity of carbon dioxide at constant pressure (J/mol/K)
!! [link](http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1#Thermo-Gas)
real(dp) function cdHeatCapacity(temp)
    real(dp), intent(in) :: temp
    real(dp) :: &
        t,&
        a=24.99735_dp,&
        b=55.18696_dp,&
        c=-33.69137_dp,&
        d=7.948387_dp,&
        e=-0.136638_dp
    t=temp/1e3
    cdHeatCapacity=a+b*t+c*t**2+d*t**3+e/t**2
end function cdHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> heat capacity of cyclo-pentane at constant pressure (J/mol/K)
!! fitted to data from [link](http://webbook.nist.gov/cgi/cbook.cgi?ID=C287923&Units=SI&Mask=1#Thermo-Gas)
real(dp) function cypHeatCapacity(temp)
    real(dp), intent(in) :: temp
    real(dp) :: &
        t,&
        a=-25.6132057_dp,&
        b=226.4176882_dp,&
        c=574.2688767_dp,&
        d=-670.5517907_dp,&
        e=0.6765321_dp
    t=temp/1e3
    cypHeatCapacity=a+b*t+c*t**2+d*t**3+e/t**2
end function cypHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> heat capacity of air at constant pressure (J/mol/K)
real(dp) function airHeatCapacity(temp)
    real(dp), intent(in) :: temp
    airHeatCapacity=0.21_dp*oxyHeatCapacity(temp)+0.79_dp*nitrHeatCapacity(temp)
end function airHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> heat capacity of nitrogen at constant pressure (J/mol/K)
!! [link](http://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas)
real(dp) function nitrHeatCapacity(temp)
    real(dp), intent(in) :: temp
    real(dp) :: &
        t,&
        a=28.98641_dp,&
        b=1.853978_dp,&
        c=-9.647459_dp,&
        d=16.63537_dp,&
        e=0.000117_dp
    t=temp/1e3
    nitrHeatCapacity=a+b*t+c*t**2+d*t**3+e/t**2
end function nitrHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> heat capacity of oxygen at constant pressure (J/mol/K)
!! [link](http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas)
real(dp) function oxyHeatCapacity(temp)
    real(dp), intent(in) :: temp
    real(dp) :: &
        t,&
        a=31.32234_dp,&
        b=-20.23531_dp,&
        c=57.86644_dp,&
        d=-36.50624_dp,&
        e=-0.007374_dp
    t=temp/1e3
    oxyHeatCapacity=a+b*t+c*t**2+d*t**3+e/t**2
end function oxyHeatCapacity
!***********************************END****************************************
end module physicalProperties
