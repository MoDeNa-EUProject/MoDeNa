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
    public polymerConductivity,gasConductivity,strutContent,get_names
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
subroutine gasConductivity(kgas,temp,xgas,gasname)
    real(dp), intent(out) :: kgas !< thermal conductivity
    real(dp), intent(in) :: temp !< temperature
    real(dp), dimension(:), intent(in) :: xgas !< molar fraction of gases
    character(len=255), dimension(:), intent(in) :: gasname !< names of gases
    integer :: i
    !modena variables
    integer(c_size_t) :: kgasTemppos
    integer(c_size_t), dimension(:), allocatable :: kgasXpos
    integer(c_int) :: ret
    type(c_ptr) :: kgasModena = c_null_ptr
    type(c_ptr) :: kgasInputs = c_null_ptr
    type(c_ptr) :: kgasOutputs = c_null_ptr
    allocate(kgasXpos(size(xgas)))
    kgasModena = modena_model_new (&
        c_char_"gasMixtureConductivity"//c_null_char)
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    kgasInputs = modena_inputs_new (kgasModena)
    kgasOutputs = modena_outputs_new (kgasModena)
    kgasTemppos = modena_model_inputs_argPos(&
        kgasModena, c_char_"T"//c_null_char)
    do i=1,size(xgas)
        kgasXpos(i) = modena_model_inputs_argPos(&
            kgasModena,&
            c_char_"x[A="//TRIM(ADJUSTL(gasname(i)))//"]"//c_null_char)
    enddo
    call modena_model_argPos_check(kgasModena)
    call modena_inputs_set(kgasInputs, kgasTemppos, temp)
    do i=1,size(xgas)
        call modena_inputs_set(kgasInputs, kgasXpos(i), xgas(i))
    enddo
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


!********************************BEGINNING*************************************
!> Get names of the species.
subroutine get_names(gasname)
    character(len=255), dimension(:), allocatable :: gasname
    integer(c_size_t) :: i
    !modena variables
    type(c_ptr) :: index_set = c_null_ptr
    type(c_ptr) :: name = c_null_ptr
    integer(c_size_t) :: itbeg,itend
    index_set = modena_index_set_new(&
        c_char_"gas_thermal_conductivity_species"//c_null_char)
    itbeg = modena_index_set_iterator_start(index_set)
    itend = modena_index_set_iterator_end(index_set)
    allocate(gasname(itend))
    do i = itbeg, itend - 1
        gasname(i+1) = modena_index_set_get_name(index_set, i)
    enddo
    call modena_index_set_destroy(index_set)
end subroutine get_names
!***********************************END****************************************


!********************************BEGINNING*************************************
function modena_index_set_get_name(indexSet, ind) result(ret)
    type(c_ptr) :: indexSet
    integer(c_size_t) :: ind
    character*255 :: ret
    type(c_ptr) :: name_ptr = c_null_ptr
    character, pointer, dimension(:) :: last_message_array
    character*255 :: last_message
    integer :: message_length, i
    name_ptr = modena_index_set_get_name_ptr(indexSet, ind)
    call C_F_POINTER(name_ptr, last_message_array, [ 255 ])
    do i=1, 255
        last_message(i:i+1) = last_message_array(i)
    enddo
    message_length = LEN_TRIM(last_message(1:INDEX(last_message, CHAR(0))))
    ret = last_message(1:message_length-1)
end function modena_index_set_get_name
!***********************************END****************************************
end module physicalProperties
