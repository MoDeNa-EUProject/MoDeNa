!> @file
!! contains definitions of modena variables
!! creates and destroys modena models
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module modenastuff
    use globals
    use iso_c_binding
    use fmodena
    implicit none
    private
    integer(c_int) :: ret
    integer(c_size_t) :: &
        viscTpos,viscXPos,&
        rhopTPos,rhopXOHPos,&
        itensTPos,&
        diffTPos(2),&
        solTPos(2),solXgasPos(2),solXmdiPos(2),solXpolyolPos(2),&
        kinInputsPos(21),kinOutputsPos(20)
    type(c_ptr) :: &
        viscModena = c_null_ptr,&
        viscInputs = c_null_ptr,&
        viscOutputs = c_null_ptr,&
        rhopModena = c_null_ptr,&
        rhopInputs = c_null_ptr,&
        rhopOutputs = c_null_ptr,&
        itensModena = c_null_ptr,&
        itensInputs = c_null_ptr,&
        itensOutputs = c_null_ptr,&
        diffModena(2) = c_null_ptr,&
        diffInputs(2) = c_null_ptr,&
        diffOutputs(2) = c_null_ptr,&
        solModena(2) = c_null_ptr,&
        solInputs(2) = c_null_ptr,&
        solOutputs(2) = c_null_ptr,&
        kinModena = c_null_ptr,&
        kinInputs = c_null_ptr,&
        kinOutputs = c_null_ptr
    public ret,&
        createModenaModels,destroyModenaModels,&
        viscTpos,viscXPos,&
        rhopTPos,rhopXOHPos,&
        itensTPos,&
        diffTPos,&
        solTPos,solXgasPos,solXmdiPos,solXpolyolPos,&
        kinInputsPos,kinOutputsPos,&
        viscModena,viscInputs,viscOutputs,&
        rhopModena,rhopInputs,rhopOutputs,&
        itensModena,itensInputs,itensOutputs,&
        diffModena,diffInputs,diffOutputs,&
        solModena,solInputs,solOutputs,&
        kinModena,kinInputs,kinOutputs
contains
!********************************BEGINNING*************************************
!> creates Modena models
subroutine createModenaModels
    integer :: i
    character(len=3) gasname(2)
    if (visc_model==3) then
        viscModena = modena_model_new (c_char_"polymerViscosity"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        viscInputs = modena_inputs_new (viscModena);
        viscOutputs = modena_outputs_new (viscModena);
        viscTpos = modena_model_inputs_argPos(viscModena, &
            c_char_"T"//c_null_char);
        viscXPos = modena_model_inputs_argPos(viscModena, &
            c_char_"X"//c_null_char);
        call modena_model_argPos_check(viscModena)
    endif
    if (rhop_model==2) then
        rhopModena = modena_model_new (&
            c_char_"density_reaction_mixture"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        rhopInputs = modena_inputs_new (rhopModena);
        rhopOutputs = modena_outputs_new (rhopModena);
        rhopTpos = modena_model_inputs_argPos(rhopModena, &
            c_char_"T"//c_null_char);
        rhopXOHPos = modena_model_inputs_argPos(rhopModena, &
            c_char_"XOH"//c_null_char);
        call modena_model_argPos_check(rhopModena)
    elseif (rhop_model==3) then
        rhopModena = modena_model_new (&
            c_char_"PolymerDensity[A=AIR,B=PU]"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        rhopInputs = modena_inputs_new (rhopModena);
        rhopOutputs = modena_outputs_new (rhopModena);
        rhopTpos = modena_model_inputs_argPos(rhopModena, &
            c_char_"T"//c_null_char);
        call modena_model_argPos_check(rhopModena)
    endif
    if (itens_model==2) then
        itensModena = modena_model_new (&
            c_char_"SurfaceTension[A=AIR,B=PU]"//c_null_char);
            ! c_char_"SurfaceTension[A=AIR,B=THF]"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        itensInputs = modena_inputs_new (itensModena);
        itensOutputs = modena_outputs_new (itensModena);
        itensTpos = modena_model_inputs_argPos(&
            itensModena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(itensModena)
    endif
    if (co2_pos==2) then
        gasname(1)="CyP"
        gasname(2)="CO2"
    elseif (co2_pos==1) then
        gasname(2)="CyP"
        gasname(1)="CO2"
    else
        stop 'CO2 position must be 1 or 2'
    endif
    do i=1,ngas
        if (diff_model(i)==2) then
            diffModena(i) = modena_model_new (&
                c_char_"diffusivityPol[A="//gasname(i)//"]"//c_null_char);
            if (modena_error_occurred()) then
                call exit(modena_error())
            endif
            diffInputs(i) = modena_inputs_new (diffModena(i));
            diffOutputs(i) = modena_outputs_new (diffModena(i));
            diffTpos(i) = modena_model_inputs_argPos(diffModena(i), &
                c_char_"T"//c_null_char);
            call modena_model_argPos_check(diffModena(i))
        endif
        if (sol_model(i)==1) then
            solModena(i) = modena_model_new (&
                c_char_"SolubilityCO2Baser"//c_null_char);
            if (modena_error_occurred()) then
                call exit(modena_error())
            endif
            solInputs(i) = modena_inputs_new (solModena(i));
            solOutputs(i) = modena_outputs_new (solModena(i));
            solTpos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"T"//c_null_char);
            call modena_model_argPos_check(solModena(i))
        endif
        if (sol_model(i)==2) then
            solModena(i) = modena_model_new (&
                c_char_"Solubility[A="//gasname(i)//",B=3]"//c_null_char);
            if (modena_error_occurred()) then
                call exit(modena_error())
            endif
            solInputs(i) = modena_inputs_new (solModena(i));
            solOutputs(i) = modena_outputs_new (solModena(i));
            solTpos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"T"//c_null_char);
            ! solXgasPos(i) = modena_model_inputs_argPos(solModena(i), &
            !     c_char_"xl1"//c_null_char);
            ! solXpolyolPos(i) = modena_model_inputs_argPos(solModena(i), &
            !     c_char_"xl2"//c_null_char);
            ! solXmdiPos(i) = modena_model_inputs_argPos(solModena(i), &
            !     c_char_"xl3"//c_null_char);
            call modena_model_argPos_check(solModena(i))
        endif
        if (sol_model(i)==3) then
            solModena(i) = modena_model_new (&
                c_char_"SolubilityPentGupta"//c_null_char);
            if (modena_error_occurred()) then
                call exit(modena_error())
            endif
            solInputs(i) = modena_inputs_new (solModena(i));
            solOutputs(i) = modena_outputs_new (solModena(i));
            solTpos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"T"//c_null_char);
            call modena_model_argPos_check(solModena(i))
        endif
        if (sol_model(i)==4) then
            solModena(i) = modena_model_new (&
                c_char_"SolubilityPentWinkler"//c_null_char);
            if (modena_error_occurred()) then
                call exit(modena_error())
            endif
            solInputs(i) = modena_inputs_new (solModena(i));
            solOutputs(i) = modena_outputs_new (solModena(i));
            solTpos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"T"//c_null_char);
            call modena_model_argPos_check(solModena(i))
        endif
        if (sol_model(i)==6) then
            solModena(i) = modena_model_new (&
                c_char_"SolubilityR11Baser"//c_null_char);
            if (modena_error_occurred()) then
                call exit(modena_error())
            endif
            solInputs(i) = modena_inputs_new (solModena(i));
            solOutputs(i) = modena_outputs_new (solModena(i));
            solTpos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"T"//c_null_char);
            call modena_model_argPos_check(solModena(i))
        endif
    enddo
    if (kin_model==4) then
        kinModena = modena_model_new (c_char_"RF-1-public"//c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        kinInputs = modena_inputs_new (kinModena);
        kinOutputs = modena_outputs_new (kinModena);
        kinInputsPos(1) = modena_model_inputs_argPos(kinModena, &
            c_char_"'kineticTime'"//c_null_char);
        kinInputsPos(2) = modena_model_inputs_argPos(kinModena, &
            c_char_"'Catalyst_1'"//c_null_char);
        kinInputsPos(3) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_A0'"//c_null_char);
        kinInputsPos(4) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_A1'"//c_null_char);
        kinInputsPos(5) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_B'"//c_null_char);
        kinInputsPos(6) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_B2'"//c_null_char);
        kinInputsPos(7) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_I0'"//c_null_char);
        kinInputsPos(8) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_I1'"//c_null_char);
        kinInputsPos(9) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_I2'"//c_null_char);
        kinInputsPos(10) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_PBA'"//c_null_char);
        kinInputsPos(11) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Breac'"//c_null_char);
        kinInputsPos(12) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Areac0'"//c_null_char);
        kinInputsPos(13) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Areac1'"//c_null_char);
        kinInputsPos(14) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Ireac0'"//c_null_char);
        kinInputsPos(15) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Ireac1'"//c_null_char);
        kinInputsPos(16) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Ireac2'"//c_null_char);
        kinInputsPos(17) = modena_model_inputs_argPos(kinModena, &
            c_char_"'Bulk'"//c_null_char);
        kinInputsPos(18) = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1'"//c_null_char);
        kinInputsPos(19) = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1_mass'"//c_null_char);
        kinInputsPos(20) = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1_temp'"//c_null_char);
        kinInputsPos(21) = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1_vol'"//c_null_char);
        kinOutputsPos(1) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_Catalyst_1"//c_null_char);
        kinOutputsPos(2) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_A0"//c_null_char);
        kinOutputsPos(3) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_A1"//c_null_char);
        kinOutputsPos(4) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_B"//c_null_char);
        kinOutputsPos(5) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_B2"//c_null_char);
        kinOutputsPos(6) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_I0"//c_null_char);
        kinOutputsPos(7) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_I1"//c_null_char);
        kinOutputsPos(8) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_I2"//c_null_char);
        kinOutputsPos(9) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_PBA"//c_null_char);
        kinOutputsPos(10) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_Breac"//c_null_char);
        kinOutputsPos(11) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_Areac0"//c_null_char);
        kinOutputsPos(12) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_Areac1"//c_null_char);
        kinOutputsPos(13) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_Ireac0"//c_null_char);
        kinOutputsPos(14) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_Ireac1"//c_null_char);
        kinOutputsPos(15) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CE_Ireac2"//c_null_char);
        kinOutputsPos(16) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_Bulk"//c_null_char);
        kinOutputsPos(17) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_R_1"//c_null_char);
        kinOutputsPos(18) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_R_1_mass"//c_null_char);
        kinOutputsPos(19) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_R_1_temp"//c_null_char);
        kinOutputsPos(20) = modena_model_outputs_argPos(kinModena, &
            c_char_"source_R_1_vol"//c_null_char);
        call modena_model_argPos_check(kinModena)
    endif
end subroutine createModenaModels
!***********************************END****************************************


!********************************BEGINNING*************************************
!> destroys Modena models
subroutine destroyModenaModels
    integer :: i
    if (visc_model==3) then
        call modena_inputs_destroy (viscInputs);
        call modena_outputs_destroy (viscOutputs);
        call modena_model_destroy (viscModena);
    endif
    if (rhop_model==2) then
        call modena_inputs_destroy (rhopInputs);
        call modena_outputs_destroy (rhopOutputs);
        call modena_model_destroy (rhopModena);
    endif
    if (itens_model==2) then
        call modena_inputs_destroy (itensInputs);
        call modena_outputs_destroy (itensOutputs);
        call modena_model_destroy (itensModena);
    endif
    do i=1,ngas
        if (diff_model(i)==2) then
            call modena_inputs_destroy (diffInputs(i));
            call modena_outputs_destroy (diffOutputs(i));
            call modena_model_destroy (diffModena(i));
        endif
    enddo
    do i=1,ngas
        if (sol_model(i)==2 .or. sol_model(i) == 6) then
            call modena_inputs_destroy (solInputs(i));
            call modena_outputs_destroy (solOutputs(i));
            call modena_model_destroy (solModena(i));
        endif
    enddo
    if (kin_model==2 .or. kin_model==4) then
        call modena_inputs_destroy (kinInputs);
        call modena_outputs_destroy (kinOutputs);
        call modena_model_destroy (kinModena);
    endif
end subroutine destroyModenaModels
!***********************************END****************************************
end module modenastuff
