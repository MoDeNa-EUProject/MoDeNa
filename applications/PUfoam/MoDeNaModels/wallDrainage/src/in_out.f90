!> @file
!! handles input and output
!! @author    Pavel Ferkl
!! @ingroup   wall_drain
module in_out
    implicit none
    private
    integer, dimension(:), allocatable :: fi
    public read_inputs,save_int_header,save_int_step,save_int_close,&
        load_bubble_growth
contains
!********************************BEGINNING*************************************
!> reads inputs from json file
subroutine read_inputs
    use globals
    use fson
    use fson_value_m, only: fson_value_get
    character(len=1024) :: strval
    type(fson_value), pointer :: json_data
    json_data => fson_parse("../inputs/wallDrainage_inputs.json")
    call fson_get(json_data, "initialConditions.centerThickness", hi)
    call fson_get(json_data, "growthRateModel", growthRateModel)
    if (growthRateModel=="constantGrowth") then
        call fson_get(json_data, "growthRate", gr)
        call fson_get(json_data, "initialConditions.domainSize", rd)
        call fson_get(json_data, "initialConditions.filmReduction", dstr)
        call fson_get(json_data, "integration.initialTime", initialTime)
        call fson_get(json_data, "integration.outerTimeSteps", its)
    elseif (growthRateModel=="fromFile") then
        continue
    else
        print*, 'unknown growth rate model'
        stop
    endif
    call fson_get(json_data, "physicalProperties.viscosityModel", viscosityModel)
    if (viscosityModel=="constant") then
        call fson_get(json_data, "physicalProperties.viscosity", mu)
    elseif (viscosityModel=="fromFile") then
        continue
    else
        print*, 'unknown viscosity model'
        stop
    endif
    call fson_get(json_data, "physicalProperties.surfaceTension", gam)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.N", ndp)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.M", mdp)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.hst", cdp)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.C", hdp)
    call fson_get(json_data, "physicalProperties.disjoiningPressure.B1", bdp)
    call fson_get(json_data, "integration.timeStep", timestep)
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
    call fson_get(json_data, "strutFilmParameter", strutFilmParameter)
    call fson_destroy(json_data)
end subroutine read_inputs
!***********************************END****************************************


!********************************BEGINNING*************************************
!> creates headers for outputs
subroutine save_int_header
    use ioutils, only: newunit
    allocate(fi(10))
    open (unit=newunit(fi(1)), file = 'filmthickness.csv')
    open (unit=newunit(fi(2)), file = 'results_1d.csv')
    write(*,'(1x,100a12)') 'time:','dr:','total: ','struts:',&
        'hmin: ','hloc: ','hcenter: ','havg: '
    write(unit=fi(2), fmt='(10000a12)') 'time','dr','np','vt','fs',&
        'hmin','hloc','hcenter','havg'
end subroutine save_int_header
!***********************************END****************************************


!********************************BEGINNING*************************************
!> saves results at current time
subroutine save_int_step(y,t)
    use constants, only: dp
    use globals, only: dr
    use phys_prop, only: volume_balance,min_film_thickness
    real(dp), intent(in) :: t
    real(dp), dimension(:), intent(in) :: y
    integer :: neq
    real(dp) :: fs,vt,hmin,hloc,havg
    neq=size(y)
    call volume_balance(y,vt,fs)
    call min_film_thickness(y,hmin,hloc,havg)
    write(*,'(100es12.3)') t,dr,vt,fs,hmin,hloc,y(1),havg
    write(fi(1),"(10000es12.4)") y(1:neq)
    write(unit=fi(2), fmt='(10000es12.4)') &
        t,dr,dble(neq),vt,fs,hmin,hloc,y(1),havg
end subroutine save_int_step
!***********************************END****************************************


!********************************BEGINNING*************************************
!> closes files
subroutine save_int_close
    close(fi(1))
    close(fi(2))
end subroutine save_int_close
!***********************************END****************************************


!********************************BEGINNING*************************************
!> loads evolution of bubble growth and
subroutine load_bubble_growth(matrix)
    use constants, only: dp
    use ioutils, only: newunit
    integer :: i,j,ios,fi
    real(dp), dimension(:,:), allocatable :: matrix
    j=0
    open(newunit(fi),file='../results/bubbleGrowth/bblgr_2_drain.out')
        do  !find number of points
            read(fi,*,iostat=ios)
            if (ios/=0) exit
            j=j+1
        enddo
        if (allocated(matrix)) deallocate(matrix)
        allocate(matrix(j-1,4))
        rewind(fi)
        read(fi,*)
        do i=1,j-1
            read(fi,*) matrix(i,:)
        enddo
    close(fi)
end subroutine load_bubble_growth
!***********************************END****************************************
end module in_out
