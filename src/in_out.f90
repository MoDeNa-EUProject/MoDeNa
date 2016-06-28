!handling of input and output
module in_out
    implicit none
    private
    integer, dimension(:), allocatable :: fi
    public read_inputs,save_int_header,save_int_step,save_int_close
contains
!********************************BEGINNING*************************************
! reads inputs from json file
subroutine read_inputs
    use globals
    use fson
    use fson_value_m, only: fson_value_get
    character(len=1024) :: strval
    type(fson_value), pointer :: json_data
    json_data => fson_parse("./inputs.json")
    call fson_get(json_data, "initialConditions.hi", hi)
    call fson_get(json_data, "initialConditions.rd", rd)
    call fson_get(json_data, "initialConditions.dstr", dstr)
    call fson_get(json_data, "growthRate", gr)
    call fson_get(json_data, "physicalProperties.viscosity", mu)
    call fson_get(json_data, "physicalProperties.surfaceTension", gam)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.N", ndp)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.M", mdp)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.hst", cdp)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.C", hdp)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.B1", bdp)
    call fson_get(json_data, "integration.initialTime", initialTime)
    call fson_get(json_data, "integration.timeStep", timestep)
    call fson_get(json_data, "integration.outerTimeSteps", its)
    call fson_get(json_data, "integration.method", strval)
    if (strval=="stiff") then
        int_method=222
    else
        print*, "unknown integration method"
        stop
    endif
    call fson_get(json_data, "integration.internalNodes", meshpoints)
    call fson_get(json_data, "integration.maxInnerTimeSteps", maxts)
    call fson_get(json_data, "integration.relativeTolerance", int_reltol)
    call fson_get(json_data, "integration.absoluteTolerance", int_abstol)
    call fson_get(json_data, "algebraicEquationSolver.tolerance", ae_tol)
    call fson_destroy(json_data)
end subroutine read_inputs
!***********************************END****************************************


!********************************BEGINNING*************************************
! creates headers for outputs
subroutine save_int_header
    use ioutils, only: newunit
    allocate(fi(10))
    open (unit=newunit(fi(1)), file = 'filmthickness.csv')
    open (unit=newunit(fi(2)), file = 'results_1d.csv')
    write(*,'(1x,100a12)') 'time:','dr:','film: ','strut: ','total: ',&
        'hmin ','hloc '
    write(unit=fi(2), fmt='(10000a12)') '#time','dr','np','vf','vs','vt',&
        'hmin','hloc'
end subroutine save_int_header
!***********************************END****************************************


!********************************BEGINNING*************************************
! saves results at current time
subroutine save_int_step(y,t)
    use constants, only: dp
    use globals, only: dr
    use phys_prop, only: volume_balance,min_film_thickness
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    integer :: neq
    real(dp) :: vf,vs,vt,hmin,hloc
    neq=size(y)
    call volume_balance(y,vf,vs,vt)
    call min_film_thickness(y,hmin,hloc)
    write(*,'(100es12.3)') t,dr,vf,vs,vt,hmin,hloc
    write(fi(1),"(10000es12.4)") y(1:neq)
    write(unit=fi(2), fmt='(10000es12.4)') t,dr,dble(neq),vf,vs,vt,hmin,hloc
end subroutine save_int_step
!***********************************END****************************************


!********************************BEGINNING*************************************
! closes files
subroutine save_int_close
    close(fi(1))
    close(fi(2))
end subroutine save_int_close
!***********************************END****************************************
end module in_out
