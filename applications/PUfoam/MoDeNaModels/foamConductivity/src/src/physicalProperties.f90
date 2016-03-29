!> @file
!! subroutines for calculation of physical properties of polymer and blowing
!! agents (using Modena calls)
!! @author    Pavel Ferkl
!! @ingroup   foam_cond
module physicalProperties
    use constants
    use fmodena
    implicit none
    private
    public polymerConductivity,gasConductivity
contains
!********************************BEGINNING*************************************
!> calculation of thermal conductivity of polymer
subroutine polymerConductivity(ksol,temp)
    real(dp), intent(out) :: ksol
    real(dp), intent(in) :: temp
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
!> calculation of thermal conductivity of gas
subroutine gasConductivity(kgas,temp,xCO2,xAir,xCyP)
    real(dp), intent(out) :: kgas
    real(dp), intent(in) :: temp,xCO2,xAir,xCyP
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
    call modena_inputs_set(kgasInputs, kgasXO2pos, xAir*0.21_dp)
    call modena_inputs_set(kgasInputs, kgasXN2pos, xAir*0.79_dp)
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
end module physicalProperties
