!handling of input and output
module in_out
    use globals
    implicit none
    private
    public read_inputs
contains
!********************************BEGINNING*************************************
! reads inputs from json file
subroutine read_inputs
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
end module in_out
