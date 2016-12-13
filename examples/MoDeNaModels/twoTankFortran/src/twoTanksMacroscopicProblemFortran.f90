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
!    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
!    details.
!
!    You should have received a copy of the GNU General Public License along
!    with Modena.  If not, see <http://www.gnu.org/licenses/>.
!
!Description
!    Solving the two tank problem the MoDeNa way.
!
!    A prototypical macros-scopic code embeds a micro-scale model (flowRate)
!    through the MoDeNa interface library.
!
!    Re-programmed to Fortran
!
!Authors
!    Henrik Rusche
!
!Contributors
!    Pavel Ferkl
program twoTanksMacroscopicProblemFortran !command line arguments not implemented yet
    use iso_c_binding
    use fmodena

    implicit none

    integer, parameter ::dp=selected_real_kind(15)

    real(dp), parameter :: D = 0.01_dp;

    real(dp) :: p0 = 3e5;
    real(dp) :: p1 = 10000;
    real(dp) :: V0 = 0.1_dp;
    real(dp) :: V1 = 1;
    real(dp) :: temp = 300;

    real(dp) :: t = 0.0_dp;
    real(dp), parameter :: deltat = 1e-3_dp;
    real(dp), parameter :: tend = 5.5_dp;

    real(dp) :: m0
    real(dp) :: m1

    real(dp) :: rho0
    real(dp) :: rho1

    real(dp) :: mdot

    integer(c_int) :: ret

    type(c_ptr) :: client = c_null_ptr !mongoc_client_t
    type(c_ptr) :: model = c_null_ptr !modena_model_t
    type(c_ptr) :: inputs = c_null_ptr !modena_inputs_t
    type(c_ptr) :: outputs = c_null_ptr !modena_outputs_t

    integer(c_size_t) :: Dpos
    integer(c_size_t) :: rho0Pos
    integer(c_size_t) :: p0Pos
    integer(c_size_t) :: p1Byp0Pos

    m0 = p0*V0/287.1_dp/temp;
    m1 = p1*V1/287.1_dp/temp;
    rho0 = m0/V0;
    rho1 = m1/V1;

!    //Instantiate a model
    model = modena_model_new (c_char_"flowRate"//c_null_char); !modena_model_t

!    //Allocate memory and fetch arg positions
    inputs = modena_inputs_new (model); !modena_inputs_t
    outputs = modena_outputs_new (model); !modena_outputs_t

    Dpos = modena_model_inputs_argPos(model, c_char_"D"//c_null_char);
    rho0Pos = modena_model_inputs_argPos(model, c_char_"rho0"//c_null_char);
    p0Pos = modena_model_inputs_argPos(model, c_char_"p0"//c_null_char);
    p1Byp0Pos = modena_model_inputs_argPos(model, c_char_"p1Byp0"//c_null_char);

!    //TODO: Add checking function that makes sure that the positions of all
!    //arguments to the model are requested
    call modena_model_argPos_check(model);

    do while(t + deltat < tend + 1e-10_dp)

        t = t + deltat;

        if(p0 > p1) then
!            //Set input vector
            call modena_inputs_set(inputs, Dpos, D);
            call modena_inputs_set(inputs, rho0Pos, rho0);
            call modena_inputs_set(inputs, p0Pos, p0);
            call modena_inputs_set(inputs, p1Byp0Pos, p1/p0);

!            // Call the model
            ret = modena_model_call (model, inputs, outputs);

!            // Terminate, if requested
            if(ret /= 0) then

                call modena_inputs_destroy (inputs);
                call modena_outputs_destroy (outputs);
                call modena_model_destroy (model);

                call exit(ret)
            endif

!            // Fetch result
            mdot = modena_outputs_get(outputs, 0_c_size_t);

            m0 = m0 - mdot*deltat;
            m1 = m1 + mdot*deltat;

        else

!            // Set input vector
            call modena_inputs_set(inputs, Dpos, D);
            call modena_inputs_set(inputs, rho0Pos, rho1);
            call modena_inputs_set(inputs, p0Pos, p1);
            call modena_inputs_set(inputs, p1Byp0Pos, p0/p1);

!            // Call the model
            ret = modena_model_call (model, inputs, outputs);

!            // Terminate, if requested
            if(ret /= 0) then

                call modena_inputs_destroy (inputs);
                call modena_outputs_destroy (outputs);
                call modena_model_destroy (model);

                call exit(ret)
            endif

!            // Fetch result
            mdot = modena_outputs_get(outputs, 0_c_size_t);

            m0 = m0 + mdot*deltat;
            m1 = m1 - mdot*deltat;
        endif

        rho0 = m0/V0;
        rho1 = m1/V1;
        p0 = m0/V0*287.1_dp*temp;
        p1 = m1/V1*287.1_dp*temp;

        write(*,*) "t = ",t," rho0 = ",rho0," p0 = ",p0," p1 = ",p1
    enddo

    call modena_inputs_destroy (inputs);
    call modena_outputs_destroy (outputs);
    call modena_model_destroy (model);

    call exit(0)
end program twoTanksMacroscopicProblemFortran
