!> @file      foamAging/src/src/physicalProperties.f90
!! @ingroup   src_mod_foamAging
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @brief     Calculates material properties of the system.
!! @details
!! Also defines all Modena variables and models.
module physicalProperties
    use constants
    use fmodena
    implicit none
    !modena variables
    integer(c_int) :: ret
    type(c_ptr) :: rhopModena = c_null_ptr
    type(c_ptr) :: rhopInputs = c_null_ptr
    type(c_ptr) :: rhopOutputs = c_null_ptr
    integer(c_size_t) :: rhopTemppos
    type(c_ptr) :: kfoamModena = c_null_ptr
    type(c_ptr) :: kfoamInputs = c_null_ptr
    type(c_ptr) :: kfoamOutputs = c_null_ptr
    integer(c_size_t) :: kfoamEpspos
    integer(c_size_t) :: kfoamDcellpos
    integer(c_size_t) :: kfoamFstrutpos
    integer(c_size_t) :: kfoamTemppos
    integer(c_size_t), allocatable :: kfoamXg(:)
    type(c_ptr) :: kgasModena = c_null_ptr
    type(c_ptr) :: kgasInputs = c_null_ptr
    type(c_ptr) :: kgasOutputs = c_null_ptr
    integer(c_size_t) :: kgasTemppos
    integer(c_size_t), allocatable :: kgasXg(:)
    type(c_ptr), dimension(:), allocatable :: sgModena
    type(c_ptr), dimension(:), allocatable :: sgInputs
    type(c_ptr), dimension(:), allocatable :: sgOutputs
    integer(c_size_t), dimension(:), allocatable :: sgTemppos
    integer(c_size_t), dimension(:), allocatable :: sgxl1pos
    integer(c_size_t), dimension(:), allocatable :: sgxl2pos
    type(c_ptr), dimension(:), allocatable :: dgModena
    type(c_ptr), dimension(:), allocatable :: dgInputs
    type(c_ptr), dimension(:), allocatable :: dgOutputs
    integer(c_size_t), dimension(:), allocatable :: dgTemppos
    type(c_ptr), dimension(:), allocatable :: kgModena
    type(c_ptr), dimension(:), allocatable :: kgInputs
    type(c_ptr), dimension(:), allocatable :: kgOutputs
    integer(c_size_t), dimension(:), allocatable :: kgTemppos
contains
!********************************BEGINNING*************************************
!> Creates Modena models.
!!
!! Names of models and inputs are hardcoded here.
subroutine createModels(ngas)
    use globals, only: solModel,diffModel,gasname
    integer :: ngas
    integer :: i
    kfoamModena = modena_model_new (c_char_"foamConductivity"//c_null_char);
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    kfoamInputs = modena_inputs_new (kfoamModena);
    kfoamOutputs = modena_outputs_new (kfoamModena);
    kfoamEpspos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"eps"//c_null_char);
    kfoamDcellpos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"dcell"//c_null_char);
    kfoamFstrutpos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"fstrut"//c_null_char);
    kfoamTemppos = modena_model_inputs_argPos(&
        kfoamModena, c_char_"T"//c_null_char);
    do i=1,ngas
        kfoamXg(i) = modena_model_inputs_argPos(kfoamModena, &
            c_char_"x["//TRIM(ADJUSTL(gasname(i)))//"]"//c_null_char);
    enddo
    call modena_model_argPos_check(kfoamModena)
    kgasModena = modena_model_new (&
        c_char_"gasMixtureConductivity"//c_null_char)
    if (modena_error_occurred()) then
        call exit(modena_error())
    endif
    kgasInputs = modena_inputs_new (kgasModena)
    kgasOutputs = modena_outputs_new (kgasModena)
    kgasTemppos = modena_model_inputs_argPos(&
        kgasModena, c_char_"T"//c_null_char)
    do i=1,ngas
        kgasXg(i) = modena_model_inputs_argPos(kgasModena, &
            c_char_"x["//TRIM(ADJUSTL(gasname(i)))//"]"//c_null_char);
    enddo
    call modena_model_argPos_check(kgasModena)
    do i=1,ngas
        kgModena(i) = modena_model_new (c_char_&
            "gas_thermal_conductivity[A="//TRIM(ADJUSTL(gasname(i)))//"]"//&
            c_null_char);
        if (modena_error_occurred()) then
            call exit(modena_error())
        endif
        kgInputs(i) = modena_inputs_new (kgModena(i));
        kgOutputs(i) = modena_outputs_new (kgModena(i));
        kgTemppos(i) = modena_model_inputs_argPos(&
            kgModena(i), c_char_"T"//c_null_char);
        call modena_model_argPos_check(kgModena(i))
        if (solModel(i)==1) then
            sgModena(i) = modena_model_new (c_char_&
                "Solubility[A="//TRIM(ADJUSTL(gasname(i)))//",B=2]"//&
                c_null_char);
            if (modena_error_occurred()) then
                call exit(modena_error())
            endif
            sgInputs(i) = modena_inputs_new (sgModena(i));
            sgOutputs(i) = modena_outputs_new (sgModena(i));
            sgTemppos(i) = modena_model_inputs_argPos(&
                sgModena(i), c_char_"T"//c_null_char);
            sgxl1pos(i) = modena_model_inputs_argPos(&
                sgModena(i), c_char_"xl1"//c_null_char);
            sgxl2pos(i) = modena_model_inputs_argPos(&
                sgModena(i), c_char_"xl2"//c_null_char);
            call modena_model_argPos_check(sgModena(i))
        endif
        if (diffModel(i)==1) then
            dgModena(i) = modena_model_new (c_char_&
                "diffusivityPol[A="//TRIM(ADJUSTL(gasname(i)))//"]"//&
                c_null_char);
            if (modena_error_occurred()) then
                call exit(modena_error())
            endif
            dgInputs(i) = modena_inputs_new (dgModena(i));
            dgOutputs(i) = modena_outputs_new (dgModena(i));
            dgTemppos(i) = modena_model_inputs_argPos(&
                dgModena(i), c_char_"T"//c_null_char);
            call modena_model_argPos_check(dgModena(i))
        endif
    enddo
end subroutine createModels
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Destroys Modena models.
!!
!! Cleans Modena models, inputs and outputs from memory.
subroutine destroyModels(ngas)
    use globals, only: solModel,diffModel
    integer, intent(in) :: ngas
    integer :: i
    call modena_inputs_destroy (kfoamInputs);
    call modena_outputs_destroy (kfoamOutputs);
    call modena_model_destroy (kfoamModena);
    call modena_inputs_destroy (kgasInputs);
    call modena_outputs_destroy (kgasOutputs);
    call modena_model_destroy (kgasModena);
    do i=1,ngas
        call modena_inputs_destroy (kgInputs(i));
        call modena_outputs_destroy (kgOutputs(i));
        call modena_model_destroy (kgModena(i));
        if (solModel(i)==1) then
            call modena_inputs_destroy (sgInputs(i));
            call modena_outputs_destroy (sgOutputs(i));
            call modena_model_destroy (sgModena(i));
        endif
        if (diffModel(i)==1) then
            call modena_inputs_destroy (dgInputs(i));
            call modena_outputs_destroy (dgOutputs(i));
            call modena_model_destroy (dgModena(i));
        endif
    enddo
end subroutine destroyModels
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculation of density of polymer.
!!
!! Modena call.
real(dp) function polymerDensity(temp)
    real(dp), intent(in) :: temp !< temperature
    call modena_inputs_set(rhopInputs, rhopTemppos, temp)
    ret = modena_model_call (rhopModena, rhopInputs, rhopOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    polymerDensity=modena_outputs_get(rhopOutputs, 0_c_size_t)
end function polymerDensity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Thermal conductivity of mixture of blowing agents.
!!
!! Modena call.
real(dp) function gasMixtureConductivity(temp,xg,ngas)
    integer, intent(in) :: ngas !< number of gases
    real(dp), intent(in) :: temp !< temperature
    real(dp), intent(in) :: xg(:) !< molar fractions  of gases
    integer :: i
    call modena_inputs_set(kgasInputs, kgasTemppos, temp)
    do i=1,ngas
        call modena_inputs_set(kgasInputs, kgasXg(i), xg(i))
    enddo
    ret = modena_model_call (kgasModena, kgasInputs, kgasOutputs)
    if(ret /= 0) then
        call exit(ret)
    endif
    gasMixtureConductivity=modena_outputs_get(kgasOutputs, 0_c_size_t)
end function gasMixtureConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Thermal conductivity of carbon dioxide.
!!
!! Modena call.
real(dp) function gasConductivity(temp,index)
    integer, intent(in) :: index !< index of gas
    real(dp), intent(in) :: temp !< temperature
    call modena_inputs_set(kgInputs(index), kgTemppos(index), temp)
    ret = modena_model_call (kgModena(index), kgInputs(index), kgOutputs(index))
    if(ret /= 0) then
        call exit(ret)
    endif
    gasConductivity=modena_outputs_get(kgOutputs(index), 0_c_size_t)
end function gasConductivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Solubility of gas.
!!
!! Modena call.
real(dp) function Solubility(temp,index)
    integer, intent(in) :: index !< index of gas
    real(dp), intent(in) :: temp !< temperature
    real(dp) :: xl1,xl2
    xl1=1.0e-3_dp
    xl2=1-xl1
    call modena_inputs_set(sgInputs(index), sgTemppos(index), temp)
    call modena_inputs_set(sgInputs(index), sgxl1pos(index), xl1)
    call modena_inputs_set(sgInputs(index), sgxl2pos(index), xl2)
    ret = modena_model_call (sgModena(index), sgInputs(index), sgOutputs(index))
    if(ret /= 0) then
        call exit(ret)
    endif
    Solubility=modena_outputs_get(sgOutputs(index), 0_c_size_t)
end function Solubility
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Diffusivity of gas.
!!
!! Modena call.
real(dp) function Diffusivity(temp,index)
    integer, intent(in) :: index !< index of gas
    real(dp), intent(in) :: temp !< temperature
    call modena_inputs_set(dgInputs(index), dgTemppos(index), temp)
    ret = modena_model_call (dgModena(index), dgInputs(index), dgOutputs(index))
    if(ret /= 0) then
        call exit(ret)
    endif
    Diffusivity=modena_outputs_get(dgOutputs(index), 0_c_size_t)
end function Diffusivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Diffusivity of gases in gas phase.
!!
!! Accoriding to Bird 1975, p.505, eq. 16.3-1.
real(dp) function gasDiffusivity(temp)
    use globals, only: pressure
    real(dp), intent(in) :: temp !< temperature
	real(dp) :: pcA, pcB, pcApcB, TcA, TcB, TcATcB
	real(dp) :: MA, MB, Mterm ,a, b, aToverTcsb
    pcA = 33.5e0_dp      !  N2
    pcB = 72.9e0_dp      ! CO2
    pcApcB = (pcA*pcB)**(1.0e0_dp/3.0e0_dp) ! CO2, N2, B-1 p. 744
    TcA = 126.2e0_dp     ! N2
    TcB = 304.2e0_dp     ! CO2
    TcATcB = (TcA*TcB)**(5.0e0_dp/12.0e0_dp)
    MA = 28.02e0_dp
    MB = 44.01e0_dp
    Mterm = dsqrt(1/MA + 1/MB)
    a = 2.7450e-4_dp ! non-polar pairs
    b = 1.823e0_dp
    aToverTcsb = a*(temp/dsqrt(TcA*TcB))**b
    ! pressure in atmospheres, cm2/s
    gasDiffusivity = (aToverTcsb*pcApcB*TcATcB*Mterm)*1.0e5_dp/pressure
    gasDiffusivity = gasDiffusivity * 1.0e-4_dp ! m2/s
end function gasDiffusivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Heat capacity of carbon dioxide at constant pressure (J/mol/K).
!!
!! [link](http://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1#Thermo-Gas)
real(dp) function cdHeatCapacity(temp)
    real(dp), intent(in) :: temp !< temperature
    real(dp) :: &
        t,&
        a=24.99735_dp,&
        b=55.18696_dp,&
        c=-33.69137_dp,&
        d=7.948387_dp,&
        e=-0.136638_dp
    t=temp/1e3
    cdHeatCapacity=a+b*t+c*t**2+d*t**3+e/t**2
end function cdHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Heat capacity of cyclo-pentane at constant pressure (J/mol/K).
!!
!! Fitted to data from [link](http://webbook.nist.gov/cgi/cbook.cgi?ID=C287923&Units=SI&Mask=1#Thermo-Gas)
real(dp) function cypHeatCapacity(temp)
    real(dp), intent(in) :: temp !< temperature
    real(dp) :: &
        t,&
        a=-25.6132057_dp,&
        b=226.4176882_dp,&
        c=574.2688767_dp,&
        d=-670.5517907_dp,&
        e=0.6765321_dp
    t=temp/1e3
    cypHeatCapacity=a+b*t+c*t**2+d*t**3+e/t**2
end function cypHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Heat capacity of air at constant pressure (J/mol/K).
!!
!! Calculated from oxygen and nitrogen.
real(dp) function airHeatCapacity(temp)
    real(dp), intent(in) :: temp !< temperature
    airHeatCapacity=0.21_dp*oxyHeatCapacity(temp)+0.79_dp*nitrHeatCapacity(temp)
end function airHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Heat capacity of nitrogen at constant pressure (J/mol/K).
!!
!! [link](http://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas)
real(dp) function nitrHeatCapacity(temp)
    real(dp), intent(in) :: temp !< temperature
    real(dp) :: &
        t,&
        a=28.98641_dp,&
        b=1.853978_dp,&
        c=-9.647459_dp,&
        d=16.63537_dp,&
        e=0.000117_dp
    t=temp/1e3
    nitrHeatCapacity=a+b*t+c*t**2+d*t**3+e/t**2
end function nitrHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Heat capacity of oxygen at constant pressure (J/mol/K).
!!
!! [link](http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas)
real(dp) function oxyHeatCapacity(temp)
    real(dp), intent(in) :: temp !< temperature
    real(dp) :: &
        t,&
        a=31.32234_dp,&
        b=-20.23531_dp,&
        c=57.86644_dp,&
        d=-36.50624_dp,&
        e=-0.007374_dp
    t=temp/1e3
    oxyHeatCapacity=a+b*t+c*t**2+d*t**3+e/t**2
end function oxyHeatCapacity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculation of strut content.
!!
!! Modena call.
subroutine strutContent(strut_content,foam_density)
    real(dp), intent(out) :: strut_content !< strut content
    real(dp), intent(in) :: foam_density !< foam density
    !modena variables
    integer(c_size_t) :: fspos
    integer(c_size_t) :: rhopos
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
!> Calculates effective diffusivity of the foam.
!!
!! [link](http://www.sciencedirect.com/science/article/pii/0017931086901481?via%3Dihub)
!! [link](http://cel.sagepub.com/cgi/doi/10.1177/0021955X02038002248)
elemental function effectiveDiffusivity(dcell,dwall,Pg,Seff,ksi) result(Deff)
    real(dp), intent(in) :: dcell !< cell size (m)
    real(dp), intent(in) :: dwall !< wall thickness (m)
    real(dp), intent(in) :: Pg !< wall permeability (mol/m/s/Pa)
    real(dp), intent(in) :: Seff !< effective solubility (mol/m3/Pa)
    real(dp), intent(in) :: ksi !< wall shape parameter
    real(dp) :: Deff !< effective diffusivity (m2/s)
    Deff = ksi*dcell/dwall*Pg/Seff
end function effectiveDiffusivity
!***********************************END****************************************


!********************************BEGINNING*************************************
!> Calculates effective solubility of the foam.
!!
!! [link](http://www.sciencedirect.com/science/article/pii/0017931086901481?via%3Dihub)
!! [link](http://cel.sagepub.com/cgi/doi/10.1177/0021955X02038002248)
elemental function effectiveSolubility(temp,Sg,eps) result(Seff)
    real(dp), intent(in) :: temp !< temperature (K)
    real(dp), intent(in) :: Sg !< solubility (mol/m3/Pa)
    real(dp), intent(in) :: eps !< foam porosity
    real(dp) :: Seff !< effective solubility (mol/m3/Pa)
    ! according to Ostrogorsky
    ! neglects dissolved gas in polymer
    Seff = 4.46e-5*1e2/1e5/temp*298.0_dp/1e-4
    ! according to Olsson
    ! including the effect of dissolved gas
    Seff = eps/Rg/temp+(1-eps)*Sg
end function effectiveSolubility
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
