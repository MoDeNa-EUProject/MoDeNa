!> @file
!! global variables (usually inputs)
!! @author    Michal Vonka
!! @author    Pavel Ferkl
!! @ingroup   foam_aging
module globals
	use constants, only: dp
	implicit none
	logical :: sheet
	integer :: nroutputs
	integer :: divwall, divcell, divsheet, ncell
	integer :: solModel(5),diffModel(4)
	real(dp) :: tend,tbeg
	real(dp) :: dcell, dwall, dfoam, dsheet
	real(dp) :: temp, temp_cond
	real(dp) :: pressure ! initial conds
	real(dp) :: fstrut,rhof,rhop,eps
	real(dp) :: Dgas
	real(dp), dimension(:), allocatable :: Sg,Dg,sheetSg,sheetDg,xg,pBg
end module globals
