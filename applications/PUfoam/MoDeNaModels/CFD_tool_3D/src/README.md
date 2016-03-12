# QmomKinetics
Personal version of QmomKinetics solver for testing and debugging

General Discreption:
--------------------
This is solver is aimed to include three kinetics related phenomenon in the foaming process.

The first equation is explaining the gelling reaction by conversion of OH group.
	ddt(rho_f*alpha_f*XOH) + div(rho_f*alpha_f*XOH*U) = rho_f*alpha_f*Q_kinOH
	Q_kinOH = ddt(XOH) = AOH*exp(-EOH/(RT))*c_OH0*(1-XOH)*(c_NCO0/c_OH0 - 2*(c_W0/c_OH0)*XW - XOH)*
			     (1/(1+Lliq*rho_PU/rho_BA))

The second equation is explaining the blowing reaction by conversion of water
	ddt(rho_f*alpha_f*XW) + div(rho_f*alpha_f*XW*U) = rho_f*alpha_f*Q_kinW
	Q_kinW = ddt(XW) = AW*exp(-EW/RT)*(1-XW)*(1/(1+Lliq*rho_PU/rho_BA))

The third equation is explaining the evaporation of the blowing agent by the load of evaporated physical blowing agent (LGas)
	ddt(rho_f*alpha_f*LGas) + div(rho_f*alpha_f*LGas*U) = rho_f*alpha_f*Q_BA
In the PU mixture we have dissolved and evaporated blowing agent and the concentration of the solved and evaporated BA
are always balanced:
		C_solved_BA = C_evaporated_BA
	Q_BA = ddt(LGas) = dLgas/dT * dT/dt + dLGas/dp * dp/dt
Thus, a correlation for the solubility of the blowing agent in PU is necessary.
	LGas(P,T) = Lliq0 - Lliq_max(P,T) : Lliq_max(P,T) < Lliq0
	LGas(P,T) = 0 			  : Lliq_max(P,T) >= Lliq0	
The initial liquid amount of BA is a known factor and we need an expression for Lliq_max(P,T).
In the Winkler PhD work the dependency on the pressure was neglected, since the pressure difference in a mold is
insignificant. The maximum solubilities were experimentally measured at different temperatures. 
Curve fitting results in:
	Lliq_max(T) = a + h*exp(-((T-T0)^2)/(2*w^2))
	a = 0.0064, h = 0.0551, T0 = 298, w = 17.8
The other approach is to use the empirical expression by Gupta for max soluble BA:
	Article: Formation of Integral Skin Polyurethane Foams, V.K, Gupta, D.V. Khakhar, Polymer and Eng Sci, 1999
	Lliq_max(P,T) = a /(exp((b-c(T-klnP))/(d-T+klnP)) - e)

constants based on Gupta paper:
a = -3.3e-4, b = 2.09e4, c = 67.5, d = 8.69e4, e = 1.07, k = 35.8

constants based on Winkler paper:
a = 3.3e-4, b = 2.09e4, c = 67.5, d = 8.69e4, e = 1.01, k = 35

Experimental data from Gupta
-----------------------------
Batch No.1
rhoPoly = 1066 kg/m3	mPoly = 119.98 g
rhoIso = 1228 kg/m3	mIso = 111.81 g
rhoCat = 1050 kg/m3	mCat = 0.08161 g
rhoBA = 626.5 kg/m3	m_n-pantane = 8.83 g
L0 = 0.037
Lliq = mBA / mPU

Mw = Eq.We * functionality (don't know functionality)

Following is how I dealt with the unknowns in the above equations:
c_OH0 ---> mPoly = 119.98 g, Eq.We = 155 g/mol, c_OH0 = 119.98 g * 1/155 mol/g = 0.7741 mol
c_NCO0 --> mIso = 111.81 g, Eq.We = 145 g/mol, c_NCO0 = 111.81 g * 1/145 mol/g = 0.7711 mol

Kind of an approximation for c_W0!
Blowing reaction in Winkler paper 2Iso + 1H2O --> Urea + CO2
				c_W0 ~ 0.5 * c_NCO0 ~ 0.5*0.7711 ~ 0.3856

Experimental data from Winkler
-------------------------------
AOH = 1
EOH = 35142.0 J/mol
R = 8.3145 J/mol K (const)

AW = 1050
EW = 27045.47 J/mol K

	
------------------------------------------------------------------------------
In this solver the moments are being calculated, too, with the following assumptions:

Assumptions:
------------
    + Computing moments at the interface
    + The foam onlye expands due to the reduction of the primary liquid density.
	This occurs bacause of (1-m1)*rho_liquid and rho_liquid is constant.
    + Thus, the above assumtion makes the foam's density independent of experimental profile of the density.
    + No mass transfer between foam and air
    + Two source terms have been considered corresponding to growth and coalescence. 
    + The kernels for growth and coalescence will be updated by the progress of the project.


