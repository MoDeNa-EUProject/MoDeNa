program workflowdummy
   use iso_c_binding
   use fmodena
!  use tfem_m
  implicit none
  real(8) :: temp, shear, conv
  real(8) :: mu, surfaceTension, Viscosity
  integer :: ret

   integer(c_size_t) :: TPos
   integer(c_size_t) :: XPos
   integer(c_size_t) :: shearPos
   type(c_ptr) :: Model, Inputs, Outputs
  print*, "Hello World"
  Model = modena_model_new (c_char_"Rheology_Arrhenius"//c_null_char);

  if ( modena_error_occurred() ) then
    call exit( modena_error() )
  end if

  print*, "Check Complete"
  Inputs = modena_inputs_new (Model);
  Outputs = modena_outputs_new (Model);
  TPos = modena_model_inputs_argPos(Model, &
             c_char_"T"//c_null_char);
  shearPos = modena_model_inputs_argPos(Model, &
             c_char_"shear"//c_null_char);
  XPos = modena_model_inputs_argPos(Model, &
             c_char_"X"//c_null_char);
!
! viscosityModel = modena_model_new (c_char_"polymerViscosity"//c_null_char);
! viscosityInputs = modena_inputs_new (viscosityModel);
! viscosityOutputs = modena_outputs_new (viscosityModel);
! viscosityTempPos = modena_model_inputs_argPos(viscosityModel, &
!            c_char_"T"//c_null_char);
! viscosityConvPos = modena_model_inputs_argPos(viscosityModel, &
!            c_char_"X"//c_null_char);

!  call modena_model_argPos_check(Model)

  print*, "Argpos"
! call modena_model_argPos_check(viscosityModel)


! open(14, file='RheologyExact.in')
! read(14, * ) temp, shear, conv
! close(14)
! shear =0.5
!  temp = 305
!  conv = 0.4
 shear =0.01
 temp = 300
 conv = 0.1
  call modena_inputs_set(Inputs, TPos, temp )
  call modena_inputs_set(Inputs, XPos, conv )
  call modena_inputs_set(Inputs, shearPos, shear )

! call modena_inputs_set( viscosityInputs, viscosityTempPos, temp )
! call modena_inputs_set( viscosityInputs, viscosityConvPos, conv )
!
  ret = modena_model_call(Model, Inputs, Outputs)
  if ( ret /= 0 ) then
    call modena_inputs_destroy( Inputs )
    call modena_outputs_destroy( Outputs )
    call modena_model_destroy( Model )
    print*,ret
    call exit(ret)
  end if

! ret = modena_model_call(viscosityModel, viscosityInputs, viscosityOutputs)
! if ( ret /= 0 ) then
!   call modena_inputs_destroy( viscosityInputs )
!   call modena_outputs_destroy( viscosityOutputs )
!   call modena_model_destroy( viscosityModel )
!   print*,ret
!   call exit(ret)
! end if

 viscosity = modena_outputs_get(Outputs, 0_c_size_t);
! polymerViscosity = modena_outputs_get(viscosityOutputs, 0_c_size_t);


! mu = 1.0
! open(15, file='RheologyExact.out')
! write(15,*) mu
! close(15)
 
 
end program workflowdummy
