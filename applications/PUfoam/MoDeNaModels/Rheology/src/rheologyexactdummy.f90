program rheologyexactdummy
!  use iso_c_binding
!  use fmodena
!  use tfem_m
  implicit none
  real(8) :: temp, shear, conv
  real(8) :: mu, surfaceTension, polymerViscosity
  integer :: ret
 
!  integer(c_size_t) :: surfaceTempPos
!  integer(c_size_t) :: viscosityTempPos, viscosityConvPos
!  type(c_ptr) :: surfaceModel, surfaceInputs, surfaceOutputs
!  type(c_ptr) :: viscosityModel, viscosityInputs, viscosityOutputs

! surfaceModel = modena_model_new (c_char_"SurfaceTension"//c_null_char);
! surfaceInputs = modena_inputs_new (surfaceModel);
! surfaceOutputs = modena_outputs_new (surfaceModel);
! surfaceTempPos = modena_model_inputs_argPos(surfaceModel, &
!            c_char_"T"//c_null_char);
! 
! viscosityModel = modena_model_new (c_char_"polymerViscosity"//c_null_char);
! viscosityInputs = modena_inputs_new (viscosityModel);
! viscosityOutputs = modena_outputs_new (viscosityModel);
! viscosityTempPos = modena_model_inputs_argPos(viscosityModel, &
!            c_char_"T"//c_null_char);
! viscosityConvPos = modena_model_inputs_argPos(viscosityModel, &
!            c_char_"X"//c_null_char);

! call modena_model_argPos_check(surfaceModel)

! call modena_model_argPos_check(viscosityModel)


  open(14, file='RheologyExact.in')
  read(14, * ) temp, shear, conv
  close(14)
 
! call modena_inputs_set( surfaceInputs, surfaceTempPos, temp )

! call modena_inputs_set( viscosityInputs, viscosityTempPos, temp )
! call modena_inputs_set( viscosityInputs, viscosityConvPos, conv )
! 
! ret = modena_model_call(surfaceModel, surfaceInputs, surfaceOutputs)
! if ( ret /= 0 ) then
!   call modena_inputs_destroy( surfaceInputs )
!   call modena_outputs_destroy( surfaceOutputs )
!   call modena_model_destroy( surfaceModel )
!   print*,ret
!   call exit(ret)
! end if

! ret = modena_model_call(viscosityModel, viscosityInputs, viscosityOutputs)
! if ( ret /= 0 ) then
!   call modena_inputs_destroy( viscosityInputs )
!   call modena_outputs_destroy( viscosityOutputs )
!   call modena_model_destroy( viscosityModel )
!   print*,ret
!   call exit(ret)
! end if

! surfaceTension = modena_outputs_get(surfaceOutputs, 0_c_size_t);
! polymerViscosity = modena_outputs_get(viscosityOutputs, 0_c_size_t);


  mu = 1.0
  open(15, file='RheologyExact.out')
  write(15,*) mu
  close(15)
 
 
end program rheologyexactdummy
