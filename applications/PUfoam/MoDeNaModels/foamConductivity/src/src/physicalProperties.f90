!> @file      foamConductivity/src/src/physicalProperties.f90
!! @ingroup   src_mod_foamConductivity
!! @author    Pavel Ferkl
!! @brief     Properties of gas and solid phase.
!! @details
!! Contains Modena models for determiantion of conductivity.
!! Also contains a Modena model for strut content.
module physicalProperties
    use constants
    use fmodena
    implicit none
    private
    public polymerConductivity,gasConductivity,strutContent
contains
!********************************BEGINNING*************************************
!> Calculation of thermal conductivity of polymer.
subroutine polymerConductivity(ksol,temp)
    real(dp), intent(out) :: ksol !< thermal conductivity
    real(dp), intent(in) :: temp !< temperature
    !modena variables
    integer(c_size_t) :: ksolTemppos

    integer(c_int) :: ret

    type(c_ptr) :: ksolModena = c_null_ptr
    type(c_ptr) :: ksolInputs = c_null_ptr
    type(c_ptr) :: ksolOutputs = c_null_ptr
    ksolModena = modena_model_new (&
        c_char_"polymer_thermal_conductivity"//c_null_char)
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    ksolInputs = modena_inputs_new (ksolModena)
    ksolOutputs = modena_outputs_new (ksolModena)
    ksolTemppos = modena_model_inputs_argPos(&
        ksolModena, c_char_"T"//c_null_char)
    call modena_model_argPos_check(ksolModena)
    call modena_inputs_set(ksolInputs, ksolTemppos, temp)
    ret = modena_model_call (ksolModena, ksolInputs, ksolOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    ksol=modena_outputs_get(ksolOutputs, 0_c_size_t)
    call modena_inputs_destroy (ksolInputs)
    call modena_outputs_destroy (ksolOutputs)
    call modena_model_destroy (ksolModena)
end subroutine polymerConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculation of thermal conductivity of gas.
subroutine gasConductivity(kgas,temp,xO2,xN2,xCO2,xCyP)
    real(dp), intent(out) :: kgas !< thermal conductivity
    real(dp), intent(in) :: temp !< temperature
    real(dp), intent(in) :: xO2 !< molar fraction of O2
    real(dp), intent(in) :: xN2 !< molar fraction of N2
    real(dp), intent(in) :: xCO2 !< molar fraction of CO2
    real(dp), intent(in) :: xCyP !< molar fraction of cyclopentane
    !modena variables
    integer(c_size_t) :: kgasTemppos
    integer(c_size_t) :: kgasXCO2pos
    integer(c_size_t) :: kgasXCyPpos
    integer(c_size_t) :: kgasXO2pos
    integer(c_size_t) :: kgasXN2pos

    integer(c_int) :: ret

    type(c_ptr) :: kgasModena = c_null_ptr
    type(c_ptr) :: kgasInputs = c_null_ptr
    type(c_ptr) :: kgasOutputs = c_null_ptr
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
    call modena_inputs_set(kgasInputs, kgasTemppos, temp)
    call modena_inputs_set(kgasInputs, kgasXCO2pos, xCO2)
    call modena_inputs_set(kgasInputs, kgasXCyPpos, xCyP)
    call modena_inputs_set(kgasInputs, kgasXO2pos, xO2)
    call modena_inputs_set(kgasInputs, kgasXN2pos, xN2)
    ret = modena_model_call (kgasModena, kgasInputs, kgasOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    kgas=modena_outputs_get(kgasOutputs, 0_c_size_t)
    call modena_inputs_destroy (kgasInputs)
    call modena_outputs_destroy (kgasOutputs)
    call modena_model_destroy (kgasModena)
end subroutine gasConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculation of strut content.
subroutine strutContent(strut_content,foam_density)
    real(dp), intent(out) :: strut_content !< strut content
    real(dp), intent(in) :: foam_density !< foam density
    !modena variables
    integer(c_size_t) :: fspos
    integer(c_size_t) :: rhopos

    integer(c_int) :: ret

    type(c_ptr) :: fsModena = c_null_ptr
    type(c_ptr) :: fsInputs = c_null_ptr
    type(c_ptr) :: fsOutputs = c_null_ptr
    fsModena = modena_model_new (c_char_"strutContent"//c_null_char)
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    fsInputs = modena_inputs_new (fsModena)
    fsOutputs = modena_outputs_new (fsModena)
    rhopos = modena_model_inputs_argPos(&
        fsModena, c_char_"rho"//c_null_char)
    call modena_model_argPos_check(fsModena)
    call modena_inputs_set(fsInputs, rhopos, foam_density)
    ret = modena_model_call (fsModena, fsInputs, fsOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    strut_content=modena_outputs_get(fsOutputs, 0_c_size_t)
    call modena_inputs_destroy (fsInputs)
    call modena_outputs_destroy (fsOutputs)
    call modena_model_destroy (fsModena)
end subroutine strutContent
!***********************************END****************************************
end module physicalProperties
