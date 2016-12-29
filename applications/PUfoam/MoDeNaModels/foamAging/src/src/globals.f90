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
	logical :: sheet !< determines if sheet is used on the outside of the foam
	integer :: nroutputs !< number of outer endme steps
	integer :: divwall !< number of grid points in wall
	integer :: divcell !< number of grid points in cell
	integer :: divsheet !< number of grid points in sheet
	integer :: ncell !< number of cells
	integer, dimension(:), allocatable :: solModel !< solubility model
	integer, dimension(:), allocatable :: diffModel !< diffusivity model
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
	real(dp), dimension(:), allocatable :: Sg !< gas solubility in polymer
	real(dp), dimension(:), allocatable :: Dg !< gas diffusiviy in polymer
	real(dp), dimension(:), allocatable :: sheetSg !< gas solubility in sheet
	real(dp), dimension(:), allocatable :: sheetDg !< gas diffusivity in sheet
	real(dp), dimension(:), allocatable :: xg !< molar fractions of gases
	real(dp), dimension(:), allocatable :: pBg !< boundary conditions
end module globals
