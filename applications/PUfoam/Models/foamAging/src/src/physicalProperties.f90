!> @file
!! subroutines for calculation of physical properties of polymer and blowing
!! agents using Modena calls
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module physicalProperties
    use constants
    use fmodena
    implicit none
    integer :: solModel(3),diffModel(3)
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
    integer(c_size_t) :: kfoamKgaspos
    integer(c_size_t) :: kfoamTemppos
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
    type(c_ptr) :: sairModena = c_null_ptr
    type(c_ptr) :: sairInputs = c_null_ptr
    type(c_ptr) :: sairOutputs = c_null_ptr
    integer(c_size_t) :: sairTemppos
    type(c_ptr) :: scypModena = c_null_ptr
    type(c_ptr) :: scypInputs = c_null_ptr
    type(c_ptr) :: scypOutputs = c_null_ptr
    integer(c_size_t) :: scypTemppos
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
!    rhopModena = modena_model_new (client, c_char_"polymerDensity"//c_null_char);
!    rhopInputs = modena_inputs_new (rhopModena);
!    rhopOutputs = modena_outputs_new (rhopModena);
!    rhopTemppos = modena_model_inputs_argPos(rhopModena, c_char_"T"//c_null_char);
!    call modena_model_argPos_check(rhopModena)
    kfoamModena = modena_model_new (c_char_"foamConductivity"//c_null_char);
    kfoamInputs = modena_inputs_new (kfoamModena);
    kfoamOutputs = modena_outputs_new (kfoamModena);
    kfoamEpspos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"eps"//c_null_char);
    kfoamDcellpos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"dcell"//c_null_char);
    kfoamFstrutpos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"fstrut"//c_null_char);
    kfoamKgaspos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"kgas"//c_null_char);
    kfoamTemppos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"T"//c_null_char);
    call modena_model_argPos_check(kfoamModena)
    kcdModena = modena_model_new (&
        c_char_"gas_thermal_conductivity[A=CO2]"//c_null_char);
    kcdInputs = modena_inputs_new (kcdModena);
    kcdOutputs = modena_outputs_new (kcdModena);
    kcdTemppos = modena_model_inputs_argPos(kcdModena, c_char_"T"//c_null_char);
    call modena_model_argPos_check(kcdModena)
    kairModena = modena_model_new (&
        c_char_"gas_thermal_conductivity[A=Air]"//c_null_char);
    kairInputs = modena_inputs_new (kairModena);
    kairOutputs = modena_outputs_new (kairModena);
    kairTemppos = modena_model_inputs_argPos(kairModena, c_char_"T"//c_null_char);
    call modena_model_argPos_check(kairModena)
    kcypModena = modena_model_new (&
        c_char_"gas_thermal_conductivity[A=CyP]"//c_null_char);
    kcypInputs = modena_inputs_new (kcypModena);
    kcypOutputs = modena_outputs_new (kcypModena);
    kcypTemppos = modena_model_inputs_argPos(kcypModena, c_char_"T"//c_null_char);
    call modena_model_argPos_check(kcypModena)
    if (solModel(1)==1) then
        sairModena = modena_model_new (&
            c_char_"solubilityPol[A=Air]"//c_null_char);
        sairInputs = modena_inputs_new (sairModena);
        sairOutputs = modena_outputs_new (sairModena);
        sairTemppos = modena_model_inputs_argPos(&
            sairModena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(sairModena)
    endif
    if (solModel(2)==1) then
        scdModena = modena_model_new (&
            c_char_"solubilityPol[A=CO2]"//c_null_char);
        scdInputs = modena_inputs_new (scdModena);
        scdOutputs = modena_outputs_new (scdModena);
        scdTemppos = modena_model_inputs_argPos(&
            scdModena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(scdModena)
    endif
    if (solModel(3)==1) then
        scypModena = modena_model_new (&
            c_char_"solubilityPol[A=CyP]"//c_null_char);
        scypInputs = modena_inputs_new (scypModena);
        scypOutputs = modena_outputs_new (scypModena);
        scypTemppos = modena_model_inputs_argPos(&
            scypModena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(scypModena)
    endif
    if (diffModel(1)==1) then
        do2Modena = modena_model_new (&
            c_char_"diffusivityPol[A=O2]"//c_null_char);
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
        dcdInputs = modena_inputs_new (dcdModena);
        dcdOutputs = modena_outputs_new (dcdModena);
        dcdTemppos = modena_model_inputs_argPos(&
            dcdModena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(dcdModena)
    endif
    if (diffModel(3)==1) then
        dcypModena = modena_model_new (&
            c_char_"diffusivityPol[A=CyP]"//c_null_char);
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
!    call modena_inputs_destroy (rhopInputs);
!    call modena_outputs_destroy (rhopOutputs);
!    call modena_model_destroy (rhopModena);
    call modena_inputs_destroy (kfoamInputs);
    call modena_outputs_destroy (kfoamOutputs);
    call modena_model_destroy (kfoamModena);
    call modena_inputs_destroy (kcdInputs);
    call modena_outputs_destroy (kcdOutputs);
    call modena_model_destroy (kcdModena);
    call modena_inputs_destroy (kairInputs);
    call modena_outputs_destroy (kairOutputs);
    call modena_model_destroy (kairModena);
    call modena_inputs_destroy (kcypInputs);
    call modena_outputs_destroy (kcypOutputs);
    call modena_model_destroy (kcypModena);
    call modena_inputs_destroy (scdInputs);
    call modena_outputs_destroy (scdOutputs);
    call modena_model_destroy (scdModena);
    call modena_inputs_destroy (sairInputs);
    call modena_outputs_destroy (sairOutputs);
    call modena_model_destroy (sairModena);
    call modena_inputs_destroy (scypInputs);
    call modena_outputs_destroy (scypOutputs);
    call modena_model_destroy (scypModena);
    call modena_inputs_destroy (dcdInputs);
    call modena_outputs_destroy (dcdOutputs);
    call modena_model_destroy (dcdModena);
    call modena_inputs_destroy (dcypInputs);
    call modena_outputs_destroy (dcypOutputs);
    call modena_model_destroy (dcypModena);
    call modena_inputs_destroy (do2Inputs);
    call modena_outputs_destroy (do2Outputs);
    call modena_model_destroy (do2Modena);
    call modena_inputs_destroy (dn2Inputs);
    call modena_outputs_destroy (dn2Outputs);
    call modena_model_destroy (dn2Modena);
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
!> solubility of carbon dioxide
real(dp) function cdSolubility(temp)
    real(dp), intent(in) :: temp
    call modena_inputs_set(scdInputs, scdTemppos, temp)
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
    call modena_inputs_set(sairInputs, sairTemppos, temp)
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
    call modena_inputs_set(scypInputs, scypTemppos, temp)
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
end module physicalProperties
