!
!
!   ooo        ooooo           oooooooooo.             ooooo      ooo
!   `88.       .888'           `888'   `Y8b            `888b.     `8'
!    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
!    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
!    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
!    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
!   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o
!
!Copyright
!    2014-2016 MoDeNa Consortium, All rights reserved.
!
!License
!    This file is part of Modena.
!
!    Modena is free software; you can redistribute it and/or modify it under
!    the terms of the GNU General Public License as published by the Free
!    Software Foundation, either version 3 of the License, or (at your option)
!    any later version.
!
!    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
!    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
!    FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!    for more details.
!
!    You should have received a copy of the GNU General Public License along
!    with Modena.  If not, see <http://www.gnu.org/licenses/>.

!>
!! Fortran bindings for MoDeNa low-level interface library
!!
!! @author Henrik Rusche
!! @author Pavel Ferkl
!! @copyright  2014-2016, MoDeNa Project. GNU Public License.

!>
!! @defgroup Fortran_interface_library
!! @{

module fmodena !command line arguments not implemented yet
use iso_c_binding

implicit none

interface
    function modena_inputs_new(model) result(inputs) bind(c)
        import
        type(c_ptr), value :: model
        type(c_ptr) :: inputs
    end function modena_inputs_new
    function modena_outputs_new(model) result(outputs) bind(c)
        import
        type(c_ptr), value :: model
        type(c_ptr) :: outputs
    end function modena_outputs_new
    function modena_model_inputs_argPos(model,name) result(argPos) bind(c,name="modena_model_inputs_argPos")
        import
        type(c_ptr), value :: model
        character(c_char) :: name(*)
        integer(c_size_t) :: argPos
    end function modena_model_inputs_argPos
    function modena_model_outputs_argPos(model,name) result(argPos) bind(c,name="modena_model_outputs_argPos")
        import
        type(c_ptr), value :: model
        character(c_char) :: name(*)
        integer(c_size_t) :: argPos
    end function modena_model_outputs_argPos
    subroutine modena_model_argPos_check(model) bind(c,name="modena_model_argPos_check")
        import
        type(c_ptr), value :: model
    end subroutine modena_model_argPos_check
    subroutine modena_inputs_set(inputs,argpos,input) bind(c)
        import
        type(c_ptr), value :: inputs
        integer(c_size_t), value :: argpos
        real(c_double), value :: input
    end subroutine modena_inputs_set
    subroutine modena_inputs_destroy(inputs) bind(c)
        import
        type(c_ptr), value :: inputs
    end subroutine modena_inputs_destroy
    subroutine modena_outputs_destroy(outputs) bind(c)
        import
        type(c_ptr), value :: outputs
    end subroutine modena_outputs_destroy
    subroutine modena_model_destroy(model) bind(c)
        import
        type(c_ptr), value :: model
    end subroutine modena_model_destroy
    function modena_outputs_get(outputs,argPos) result(output) bind(c)
        import
        type(c_ptr), value :: outputs
        integer(c_size_t), value :: argPos
        real(c_double) :: output
    end function modena_outputs_get
    function modena_model_new(model_name) result(model) bind(c)
        import
        character(c_char) :: model_name(*)
        type(c_ptr) :: model
    end function modena_model_new
    function modena_model_call(model,inputs,outputs) result(ret) bind(c)
        import
        type(c_ptr), value :: model
        type(c_ptr), value :: inputs
        type(c_ptr), value :: outputs
        integer(c_int) :: ret
    end function modena_model_call
    function modena_error_occurred() result(output) bind(c)
        import
        logical(c_bool) :: output
    end function modena_error_occurred
    function modena_error() result(output) bind(c)
        import
        integer(c_int) :: output
    end function modena_error
    function modena_index_set_new(indexSetId) result(indexSet) bind(c)
        import
        character(c_char) :: indexSetId(*)
        type(c_ptr) :: indexSet
    end function modena_index_set_new
    function modena_index_set_iterator_start(indexSet) result(ret) bind(c)
        import
        type(c_ptr), value :: indexSet
        integer(c_size_t) :: ret
    end function modena_index_set_iterator_start
    function modena_index_set_iterator_end(indexSet) result(ret) bind(c)
        import
        type(c_ptr), value :: indexSet
        integer(c_size_t) :: ret
    end function modena_index_set_iterator_end
    function modena_index_set_get_name_ptr(indexSet, ind) result(ret) bind(c, name='modena_index_set_get_name')
        import
        type(c_ptr), value :: indexSet
        integer(c_size_t), value :: ind
        type(c_ptr) :: ret
    end function modena_index_set_get_name_ptr
    subroutine modena_index_set_destroy(indexSet) bind(c)
        import
        type(c_ptr), value :: indexSet
    end subroutine modena_index_set_destroy
end interface

contains
!********************************BEGINNING*************************************
!! put this function directly to your code, I can't link it for some reason
! function fmodena_index_set_get_name(indexSet, ind) result(ret)
!     use  iso_c_binding
!     type(c_ptr) :: indexSet
!     integer(c_size_t) :: ind
!     character*255 :: ret
!     type(c_ptr) :: name_ptr = c_null_ptr
!     character, pointer, dimension(:) :: last_message_array
!     character*255 :: last_message
!     integer :: message_length, i
!     name_ptr = modena_index_set_get_name_ptr(indexSet, ind)
!     call C_F_POINTER(name_ptr, last_message_array, [ 255 ])
!     do i=1, 255
!         last_message(i:i+1) = last_message_array(i)
!     enddo
!     message_length = LEN_TRIM(last_message(1:INDEX(last_message, CHAR(0))))
!     ret = last_message(1:message_length-1)
! end function fmodena_index_set_get_name
!***********************************END****************************************
end module fmodena

!>
!!@}
