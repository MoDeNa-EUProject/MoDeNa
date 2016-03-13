!> @file
!! contains definitions of modena variables
!! creates and destroys modena models
!! @author    Pavel Ferkl
!! @ingroup   bblgr
module modenastuff
    use in_out
    use iso_c_binding
    use fmodena
    implicit none
    integer(c_int) :: ret
    integer(c_size_t) :: viscTpos
    integer(c_size_t) :: viscXPos
    integer(c_size_t) :: rhopTPos
    integer(c_size_t) :: rhopXOHPos
    integer(c_size_t) :: itensTPos
    integer(c_size_t) :: diffTPos(2)
    integer(c_size_t) :: solTPos(2)
    integer(c_size_t) :: solXgasPos(2)
    integer(c_size_t) :: solXmdiPos(2)
    integer(c_size_t) :: solXpolyolPos(2)
    integer(c_size_t) :: kinNCOPos
    integer(c_size_t) :: kinOHPos
    integer(c_size_t) :: kinH2OPos
    integer(c_size_t) :: kinCO2Pos
    integer(c_size_t) :: kinPentanePos
    integer(c_size_t) :: kinPolymerPos
    integer(c_size_t) :: kinPolymerBlowPos
    integer(c_size_t) :: kinUreaPos
    integer(c_size_t) :: kinR1Pos
    integer(c_size_t) :: kinRmassPos
    integer(c_size_t) :: kinRvolPos
    integer(c_size_t) :: kinRtempPos
    integer(c_size_t) :: kinSourceNCOPos
    integer(c_size_t) :: kinSourceOHPos
    integer(c_size_t) :: kinSourceH2OPos
    integer(c_size_t) :: kinSourceCO2Pos
    integer(c_size_t) :: kinSourcePentanePos
    integer(c_size_t) :: kinSourcePolymerPos
    integer(c_size_t) :: kinSourcePolymerBlowPos
    integer(c_size_t) :: kinSourceUreaPos
    integer(c_size_t) :: kinSourceR1Pos
    integer(c_size_t) :: kinSourceRmassPos
    integer(c_size_t) :: kinSourceRvolPos
    integer(c_size_t) :: kinSourceRtempPos
    integer(c_size_t) :: kinInputsPos(20)
    integer(c_size_t) :: kinOutputsPos(20)
    type(c_ptr) :: viscModena = c_null_ptr
    type(c_ptr) :: viscInputs = c_null_ptr
    type(c_ptr) :: viscOutputs = c_null_ptr
    type(c_ptr) :: rhopModena = c_null_ptr
    type(c_ptr) :: rhopInputs = c_null_ptr
    type(c_ptr) :: rhopOutputs = c_null_ptr
    type(c_ptr) :: itensModena = c_null_ptr
    type(c_ptr) :: itensInputs = c_null_ptr
    type(c_ptr) :: itensOutputs = c_null_ptr
    type(c_ptr) :: diffModena(2) = c_null_ptr
    type(c_ptr) :: diffInputs(2) = c_null_ptr
    type(c_ptr) :: diffOutputs(2) = c_null_ptr
    type(c_ptr) :: solModena(2) = c_null_ptr
    type(c_ptr) :: solInputs(2) = c_null_ptr
    type(c_ptr) :: solOutputs(2) = c_null_ptr
    type(c_ptr) :: kinModena = c_null_ptr
    type(c_ptr) :: kinInputs = c_null_ptr
    type(c_ptr) :: kinOutputs = c_null_ptr
    character(len=3) gasname(2)
contains
!********************************BEGINNING*************************************
!> creates Modena models
subroutine createModenaModels
    integer :: i
    if (visc_model==3) then
        viscModena = modena_model_new (c_char_"polymerViscosity"//c_null_char);
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
        rhopInputs = modena_inputs_new (rhopModena);
        rhopOutputs = modena_outputs_new (rhopModena);
        rhopTpos = modena_model_inputs_argPos(rhopModena, &
            c_char_"T"//c_null_char);
        rhopXOHPos = modena_model_inputs_argPos(rhopModena, &
            c_char_"XOH"//c_null_char);
        call modena_model_argPos_check(rhopModena)
    endif
    if (itens_model==2) then
        itensModena = modena_model_new (&
            c_char_"SurfaceTension[A=CO2,B=PU]"//c_null_char);
        itensInputs = modena_inputs_new (itensModena);
        itensOutputs = modena_outputs_new (itensModena);
        itensTpos = modena_model_inputs_argPos(&
            itensModena, c_char_"T"//c_null_char);
        call modena_model_argPos_check(itensModena)
    endif
    if (diff_model(1)==2) diffModena(1) = modena_model_new (&
        c_char_"gas_diffusivity[A=CyP]"//c_null_char);
    if (diff_model(2)==2) diffModena(2) = modena_model_new (&
        c_char_"gas_diffusivity[A=CO2]"//c_null_char);
    do i=1,ngas
        if (diff_model(i)==2) then
            diffInputs(i) = modena_inputs_new (diffModena(i));
            diffOutputs(i) = modena_outputs_new (diffModena(i));
            diffTpos(i) = modena_model_inputs_argPos(diffModena(i), &
                c_char_"T"//c_null_char);
            call modena_model_argPos_check(diffModena(i))
        endif
    enddo
    gasname(1)="CyP"
    gasname(2)="CO2"
    if (sol_model(1)==2) solModena(1) = modena_model_new (&
        c_char_"Solubility[A=CyP,B=3]"//c_null_char);
    if (sol_model(2)==2) solModena(1) = modena_model_new (&
        c_char_"Solubility[A=CO2,B=3]"//c_null_char);
    do i=1,ngas
        if (sol_model(i)==2) then
            solInputs(i) = modena_inputs_new (solModena(i));
            solOutputs(i) = modena_outputs_new (solModena(i));
            solTpos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"T"//c_null_char);
            solXgasPos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"xl["//gasname//"]"//c_null_char);
            solXmdiPos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"xl[mdi]"//c_null_char);
            solXpolyolPos(i) = modena_model_inputs_argPos(solModena(i), &
                c_char_"xl[polyol]"//c_null_char);
            call modena_model_argPos_check(solModena(i))
        endif
    enddo
    if (kin_model==2) then
        kinModena = modena_model_new (c_char_"simpleKinetics"//c_null_char);
        kinInputs = modena_inputs_new (kinModena);
        kinOutputs = modena_outputs_new (kinModena);

        kinNCOPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'EG_NCO'"//c_null_char);
        kinOHPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'EG_OH'"//c_null_char);
        kinH2OPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'H2O'"//c_null_char);
        kinCO2Pos = modena_model_inputs_argPos(kinModena, &
            c_char_"'CO2'"//c_null_char);
        kinPentanePos = modena_model_inputs_argPos(kinModena, &
            c_char_"'PENTANE'"//c_null_char);
        kinPolymerPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'POLYMER'"//c_null_char);
        kinPolymerBlowPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'POLMERBLOW'"//c_null_char);
        kinUreaPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'UREA'"//c_null_char);
        kinR1Pos = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1'"//c_null_char);
        kinRmassPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1_mass'"//c_null_char);
        kinRvolPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1_vol'"//c_null_char);
        kinRtempPos = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1_temp'"//c_null_char);

        kinSourceNCOPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_EG_NCO"//c_null_char);
        kinSourceOHPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_EG_OH"//c_null_char);
        kinSourceH2OPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_H2O"//c_null_char);
        kinSourceCO2Pos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_CO2"//c_null_char);
        kinSourcePentanePos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_PENTANE"//c_null_char);
        kinSourcePolymerPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_POLYMER"//c_null_char);
        kinSourcePolymerBlowPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_POLMERBLOW"//c_null_char);
        kinSourceUreaPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_UREA"//c_null_char);
        kinSourceR1Pos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_R_1"//c_null_char);
        kinSourceRmassPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_R_1_mass"//c_null_char);
        kinSourceRvolPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_R_1_vol"//c_null_char);
        kinSourceRtempPos = modena_model_outputs_argPos(kinModena, &
            c_char_"source_R_1_temp"//c_null_char);
        call modena_model_argPos_check(kinModena)
    elseif (kin_model==4) then
        kinModena = modena_model_new (c_char_"RF-1-public"//c_null_char);
        kinInputs = modena_inputs_new (kinModena);
        kinOutputs = modena_outputs_new (kinModena);
        kinInputsPos(1) = modena_model_inputs_argPos(kinModena, &
            c_char_"'Catalyst_1'"//c_null_char);
        kinInputsPos(2) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_A0'"//c_null_char);
        kinInputsPos(3) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_A1'"//c_null_char);
        kinInputsPos(4) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_B'"//c_null_char);
        kinInputsPos(5) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_B2'"//c_null_char);
        kinInputsPos(6) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_I0'"//c_null_char);
        kinInputsPos(7) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_I1'"//c_null_char);
        kinInputsPos(8) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_I2'"//c_null_char);
        kinInputsPos(9) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_PBA'"//c_null_char);
        kinInputsPos(10) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Breac'"//c_null_char);
        kinInputsPos(11) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Areac0'"//c_null_char);
        kinInputsPos(12) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Areac1'"//c_null_char);
        kinInputsPos(13) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Ireac0'"//c_null_char);
        kinInputsPos(14) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Ireac1'"//c_null_char);
        kinInputsPos(15) = modena_model_inputs_argPos(kinModena, &
            c_char_"'CE_Ireac2'"//c_null_char);
        kinInputsPos(16) = modena_model_inputs_argPos(kinModena, &
            c_char_"'Bulk'"//c_null_char);
        kinInputsPos(17) = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1'"//c_null_char);
        kinInputsPos(18) = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1_mass'"//c_null_char);
        kinInputsPos(19) = modena_model_inputs_argPos(kinModena, &
            c_char_"'R_1_temp'"//c_null_char);
        kinInputsPos(20) = modena_model_inputs_argPos(kinModena, &
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
        if (sol_model(i)==2) then
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
