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
	real(dp) :: pBCair, pBCCO2, pBCcyp 	! boundary conds
	real(dp) :: pICair, pICCO2, pICcyp 	! initial conds
	real(dp) :: fstrut,rhof,rhop,eps
	real(dp) :: DCO2, Dcyp, Dair, Dgas
	real(dp) :: SCO2, Scyp, Sair
	real(dp) :: sheetDCO2, sheetDcyp, sheetDair
	real(dp) :: sheetSCO2, sheetScyp, sheetSair
end module globals
