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
	real(dp) :: pBCair, pBCCO2, pBCcyp, pBCO2, pBCN2 	! boundary conds
	real(dp) :: xAir, xCO2, xCyP, xO2, xN2 	! initial conds
	real(dp) :: fstrut,rhof,rhop,eps
	real(dp) :: DCO2, Dcyp, Dair, Dgas, DO2, DN2
	real(dp) :: SCO2, Scyp, Sair, SO2, SN2
	real(dp) :: sheetDCO2, sheetDcyp, sheetDair, sheetDO2, sheetDN2
	real(dp) :: sheetSCO2, sheetScyp, sheetSair, sheetSO2, sheetSN2
end module globals
