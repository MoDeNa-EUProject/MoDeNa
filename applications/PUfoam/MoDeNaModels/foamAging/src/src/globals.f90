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
	integer :: solModel(3),diffModel(3)
	real(dp) :: tend,tbeg
	real(dp) :: dcell, dwall, dfoam, dsheet
	real(dp) :: temp, temp_cond
	real(dp) :: pressure ! initial conds
	real(dp) :: pBCair, pBCCO2, pBCpent 	! boundary conds
	real(dp) :: pICair, pICCO2, pICpent 	! initial conds
	real(dp) :: fstrut,rhof,rhop
	real(dp) :: DCO2, Dpent, Dair, Dgas
	real(dp) :: SCO2, Spent, Sair
	real(dp) :: sheetDCO2, sheetDpent, sheetDair
	real(dp) :: sheetSCO2, sheetSpent, sheetSair
end module globals
