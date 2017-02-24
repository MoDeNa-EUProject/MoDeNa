program workflowdummy
   use iso_c_binding
   use fmodena
!  use tfem_m
  implicit none
  real(8) :: temp, shear, conv, m0, m1
  real(8) :: mu, surfaceTension, Viscosity
  integer :: ret,i

   integer(c_size_t) :: TPos
   integer(c_size_t) :: XPos
   integer(c_size_t) :: shearPos
   integer(c_size_t) :: m0Pos
   integer(c_size_t) :: m1Pos
   type(c_ptr) :: Model, Inputs, Outputs
  print*, "Hello World"
  Model = modena_model_new (c_char_"Rheology_Arrhenius"//c_null_char);

  print*, "Model loaded from db"
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
  m0Pos = modena_model_inputs_argPos(Model, &
             c_char_"m0"//c_null_char);
  m1Pos = modena_model_inputs_argPos(Model, &
             c_char_"m1"//c_null_char);
!
! viscosityModel = modena_model_new (c_char_"polymerViscosity"//c_null_char);
! viscosityInputs = modena_inputs_new (viscosityModel);
! viscosityOutputs = modena_outputs_new (viscosityModel);
! viscosityTempPos = modena_model_inputs_argPos(viscosityModel, &
!            c_char_"T"//c_null_char);
! viscosityConvPos = modena_model_inputs_argPos(viscosityModel, &
!            c_char_"X"//c_null_char);

 ! call modena_model_argPos_check(Model)

  print*, "Argpos"
! call modena_model_argPos_check(viscosityModel)


! open(14, file='RheologyExact.in')
! read(14, * ) temp, shear, conv
! close(14)
! shear =0.5
!  temp = 305
!  conv = 0.4
 shear =0.02
 temp = 320
 conv = 0.8
 m0=1e10
 m1=0.2
  call modena_inputs_set(Inputs, TPos, temp )
  call modena_inputs_set(Inputs, XPos, conv )
  call modena_inputs_set(Inputs, shearPos, shear )
  call modena_inputs_set(Inputs, m0Pos, m0 )
  call modena_inputs_set(Inputs, m1Pos, m1 )

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
write(unit=*, fmt=*) viscosity

open(11,file='visc.out')
conv=0
do i=1,11
    call modena_inputs_set(Inputs, XPos, conv )
    ret = modena_model_call(Model, Inputs, Outputs)
    viscosity = modena_outputs_get(Outputs, 0_c_size_t);
    write(11,*) conv,viscosity
    conv=conv+0.4/10
enddo
close(11)
! mu = 1.0
! open(15, file='RheologyExact.out')
! write(15,*) mu
! close(15)


end program workflowdummy
