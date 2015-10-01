!subroutines for calculation of physical properties of polymer and blowing agents (using Modena calls)
!author: pavel.ferkl@vscht.cz
module physicalProperties
    use constants
    use fmodena
    implicit none
    private
    public polymerConductivity
contains
!********************************BEGINNING*************************************
!calculation of thermal conductivity of polymer
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
end module physicalProperties
