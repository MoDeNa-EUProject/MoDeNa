!> @file      foamAging/src/src/globals.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_foamAging
!! @brief     Global variables.
!! @details
!! Stores variables shared by many modules.
module globals
	use constants, only: dp
	implicit none
	logical :: sheet !< determines if sheet is used on the outside of the foam
	character(len=80) :: modelType !< heterogeneous or homogeneous
	character(len=80) :: progressTime !< linear or logarithmic
	!> names of the gases from the index set
    character(len=255), dimension(:), allocatable :: gasname(:)
	integer :: nroutputs !< number of outer endme steps
	integer :: outputsPerOrder !< outputs per decadic order in logarithmic scale
	integer :: numberOfOrders !< number of decadic orders in logarithmic output
	integer :: divwall !< number of grid points in wall
	integer :: divcell !< number of grid points in cell
	integer :: divfoam !< number of grid points in foam
	integer :: divsheet !< number of grid points in sheet
	integer :: ncell !< number of cells
	integer :: ngas !< number of gases
	integer :: nfv !< number of finite volumes
	integer, dimension(:), allocatable :: solModel !< solubility model
	integer, dimension(:), allocatable :: diffModel !< diffusivity model
	integer, dimension(:), allocatable :: mor !< 1=cell,2=wall,3=sheet
	real(dp) :: tend !< time at the beginning
	real(dp) :: tbeg !< time at the end
	real(dp) :: dcell !< cell size
	real(dp) :: dwall !< wall thickness
	real(dp) :: dfoam !< foam thickness
	real(dp) :: dsheet !< sheet thickness
	real(dp) :: temp !< temperature of aging
	real(dp) :: temp_cond !< temperature of conductivity measurements
	real(dp) :: pressure !< initial pressure
	real(dp) :: fstrut !< strut content
	real(dp) :: rhof !< foam density
	real(dp) :: rhop !< polymer density
	real(dp) :: eps !< foam porosity
	real(dp) :: Dgas !< gas diffusivity in gasphase
	real(dp) :: ksi !< wall shape parameter
	real(dp), dimension(:), allocatable :: dz !< size of finite volume
	real(dp), dimension(:), allocatable :: dif !< diffusivity
	real(dp), dimension(:), allocatable :: sol !< solubility
	real(dp), dimension(:), allocatable :: bc !< boundary condition
	real(dp), dimension(:), allocatable :: Sg !< gas solubility in polymer
	real(dp), dimension(:), allocatable :: Dg !< gas diffusivity in polymer
	real(dp), dimension(:), allocatable :: Pg !< gas permeability in polymer
	real(dp), dimension(:), allocatable :: Deff !< gas diffusivity in foam
	real(dp), dimension(:), allocatable :: Seff !< gas solubility in foam
	real(dp), dimension(:), allocatable :: sheetSg !< gas solubility in sheet
	real(dp), dimension(:), allocatable :: sheetDg !< gas diffusivity in sheet
	real(dp), dimension(:), allocatable :: xg !< molar fractions of gases
	real(dp), dimension(:), allocatable :: pBg !< boundary conditions
	real(dp), dimension(:), allocatable :: Mg !< molar mass
end module globals
