!******************************************************BEGINNING***************************************************************
!contains model subroutines
module model
!*****************************************************DECLARATION**************************************************************
    use ioutils
    use fmodena
    implicit none

	integer, parameter:: dp=selected_real_kind(15)

    !time integration variables for lsode
    integer :: IOUT
    real(dp), dimension(:), allocatable :: RWORK,Y
    integer, dimension(:), allocatable :: IWORK
    real(dp) :: JAC,TOUT,RTOL,ATOL,T
    integer :: IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NNZ, LENRAT, MF, NEQ

    !modena variables
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

    integer(c_int) :: ret

    type(c_ptr) :: kinModena = c_null_ptr
    type(c_ptr) :: kinInputs = c_null_ptr
    type(c_ptr) :: kinOutputs = c_null_ptr

    real(dp) :: R1,Rmass,Rvol,Rtemp
    real(dp) :: timestep,TEND
    integer :: its,maxts,pt,&
        kin_model=4 !2=simpleKinetics,4=RF-1-private
    integer, dimension(:), allocatable :: kineq

    logical, parameter :: logtimesteps=.true. !save output in logarithmic or linear time intervals
!*********************************************************BODY*****************************************************************
contains
!****************************BEGINNING*******************************
!main model
SUBROUTINE  kinetics
!***************************DECLARATION******************************
!******************************BODY**********************************
    call preproc
    call integ
END subroutine kinetics
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
!source terms
SUBROUTINE  FEX (NEQ, T, Y, YDOT)
!***************************DECLARATION******************************
    INTEGER :: NEQ
    real(dp) ::  T, Y(NEQ), YDOT(NEQ)
    integer :: i
!******************************BODY**********************************
    if (kin_model==2) then
        ! call modena_inputs_set(kinInputs, kinNCOPos, Y(kineq(1)));
    	! call modena_inputs_set(kinInputs, kinOHPos, Y(kineq(2)));
    	! call modena_inputs_set(kinInputs, kinH2OPos, Y(kineq(3)));
        ! call modena_inputs_set(kinInputs, kinCO2Pos, Y(kineq(4)));
        ! call modena_inputs_set(kinInputs, kinPentanePos, Y(kineq(5)));
        ! call modena_inputs_set(kinInputs, kinPolymerPos, Y(kineq(6)));
        ! call modena_inputs_set(kinInputs, kinPolymerBlowPos, Y(kineq(7)));
        ! call modena_inputs_set(kinInputs, kinUreaPos, Y(kineq(8)));
        ! call modena_inputs_set(kinInputs, kinR1Pos, Y(kineq(9)));
        ! call modena_inputs_set(kinInputs, kinRmassPos, Y(kineq(10)));
        ! call modena_inputs_set(kinInputs, kinRvolPos, Y(kineq(11)));
        ! call modena_inputs_set(kinInputs, kinRtempPos, Y(kineq(12)));
    elseif (kin_model==4) then
        do i=1,20
            call modena_inputs_set(kinInputs, kinInputsPos(i), Y(i))
        enddo
    endif
    ret = modena_model_call (kinModena, kinInputs, kinOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    if (kin_model==2) then
        YDOT(1) = modena_outputs_get(kinOutputs, kinSourceNCOPos);
        YDOT(2) = modena_outputs_get(kinOutputs, kinSourceOHPos);
        YDOT(3) = modena_outputs_get(kinOutputs, kinSourceH2OPos);
        YDOT(4) = modena_outputs_get(kinOutputs, kinSourceCO2Pos);
        YDOT(5) = modena_outputs_get(kinOutputs, kinSourcePentanePos);
        YDOT(6) = modena_outputs_get(kinOutputs, kinSourcePolymerPos);
        YDOT(7) = modena_outputs_get(kinOutputs, &
            kinSourcePolymerBlowPos);
        YDOT(8) = modena_outputs_get(kinOutputs, kinSourceUreaPos);
        YDOT(9) = modena_outputs_get(kinOutputs, kinSourceR1Pos);
        YDOT(10) = modena_outputs_get(kinOutputs, kinSourceRmassPos);
        YDOT(11) = modena_outputs_get(kinOutputs, kinSourceRvolPos);
        YDOT(12) = modena_outputs_get(kinOutputs, kinSourceRtempPos);
    elseif (kin_model==4) then
        YDOT=0
        do i=1,20
            YDOT(i) = modena_outputs_get(kinOutputs, kinOutputsPos(i))
            ! if (YDOT(2)>1e-12) then
            !     write(*,*) i,YDOT(2)
            !     stop
            ! endif
        enddo
    endif
END subroutine FEX
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
!creates Modena models
subroutine createModenaModels
!***************************DECLARATION******************************
!******************************BODY**********************************
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
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
!prepares integration
subroutine preproc
!***************************DECLARATION******************************
    integer :: i,j
!******************************BODY**********************************
    write(*,*) 'preparing simulation...'
    !determine number of equations and their indexes
    if (kin_model==2) then
        NEQ=12
        allocate(Y(NEQ))
        Y=0
        Y(1)=5
        Y(2)=5
        Y(3)=0.2_dp
        R1=1
        Rmass=1e0_dp
        Rvol=1e0_dp
        Rtemp=300e0_dp
        Y(9)=R1
        Y(10)=Rmass
        Y(11)=Rvol
        Y(12)=Rtemp
    elseif (kin_model==4) then
        NEQ=20
        allocate(Y(NEQ),kineq(NEQ))
        do i=1,neq
            kineq(i)=i
        enddo
        Y(kineq(1)) = 6.73000e-02_dp
        Y(kineq(2)) = 1.92250e+00_dp
        Y(kineq(3)) = 2.26920e+00_dp
        Y(kineq(4)) = 0.00000e+00_dp
        Y(kineq(5)) = 5.46200e-01_dp
        Y(kineq(6)) = 2.19790e+00_dp
        Y(kineq(7)) = 1.64000e+00_dp
        Y(kineq(8)) = 1.71030e+00_dp
        Y(kineq(9)) = 0.00000e+00_dp
        Y(kineq(10)) = 0.00000e+00_dp
        Y(kineq(11)) = 0.00000e+00_dp
        Y(kineq(12)) = 0.00000e+00_dp
        Y(kineq(13)) = 0.00000e+00_dp
        Y(kineq(14)) = 0.00000e+00_dp
        Y(kineq(15)) = 0.00000e+00_dp
        Y(kineq(16)) = 4.45849e+00_dp
        Y(kineq(17)) = 0.00000e+00_dp
        Y(kineq(18)) = 1.00000e+00_dp
        Y(kineq(19)) = 60!2.27000e+01_dp
        Y(kineq(20)) = 1e0_dp!8.46382e-01_dp
    endif
    !set initial values
    call createModenaModels
    T=0
    TEND=3e2_dp
    if (logtimesteps) then
        pt=16
        its=4*pt+1
    else
        timestep=(TEND-T)/its
        its=100
    endif
    !choose and set integrator
    MF=222
    ATOL=1e-8_dp
    RTOL=1e-8_dp
    maxts=500
    select case(MF)
    case(10)
        allocate(RWORK(20+16*NEQ),IWORK(30))
    case(222)
        NNZ=NEQ**2 !I really don't know, smaller numbers make problems for low p
        LENRAT=2 !depends on dp
        allocate(RWORK(int(20+(2+1._dp/LENRAT)*NNZ+(11+9._dp/LENRAT)*NEQ)),IWORK(30))
    case default
        stop 'unknown MF'
    end select
    ITASK = 1
    ISTATE = 1
    IOPT = 0
    LRW = size(RWORK)
    LIW = size(IWORK)
    IWORK(6)=maxts
    TOUT =T
    ITOL = 1 !don't change, or you must declare ATOL as ATOL(NEQ)
    write(*,*) 'done: simulation prepared'
    write(*,*)
end subroutine preproc
!********************************************************************
!*******************************END**********************************


!****************************BEGINNING*******************************
!performs integration
subroutine integ()
!***************************DECLARATION******************************
	integer :: fi
!******************************BODY**********************************
    write(*,*) 'integrating...'
    if (kin_model==2) then
        open(newunit(fi),file='conc.out')
        write(fi,'(1000A23)') '#time','EG_NCO','EG_OH','H2O','CO2','PENTANE','POLYMER','POLMERBLOW','UREA','R_1','R_1_mass','R_1_vol','R_1_temp'
    elseif(kin_model==4) then
        open(newunit(fi),file='kinetics.out')
        write(fi,'(1000A23)') "time","Catalyst_1","CE_A0","CE_A1","CE_B","CE_B2","CE_I0","CE_I1","CE_I2","CE_PBA","CE_Breac","CE_Areac0","CE_Areac1","CE_Ireac0","CE_Ireac1","CE_Ireac2","Bulk","R_1","R_1_mass","R_1_temp","R_1_vol"
    endif
    DO IOUT = 1,its
        if (logtimesteps) then
            TOUT = 10**((IOUT-its+pt*log10(TEND))/pt)
        else
            TOUT = TOUT+timestep
        endif
        call DLSODES (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
        write(fi,'(100es23.15)') TOUT,Y
    END DO
    close(fi)
    write(*,*) 'done: integration'

    call modena_inputs_destroy (kinInputs);
    call modena_outputs_destroy (kinOutputs);
    call modena_model_destroy (kinModena);
    call exit(0)
end subroutine integ
!********************************************************************
!*******************************END**********************************
end module model
!**********************************************************END*****************************************************************
