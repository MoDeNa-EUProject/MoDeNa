!> @file      bubbleGrowth/src/src/globals.f90
!! @author    Pavel Ferkl
!! @ingroup   src_mod_bubbleGrowth
!! @brief     Stores global variables.
!! @details
!! Defines variables, which are used throughout the program by many functions.
!! Usually they are not passed as arguments.
module globals
    use constants
    implicit none
    character(len=99) :: &
        fileplacein,& !< location of input files
        fileplaceout,& !< location of output files
        inputs,& !< input file
        outputs_1d,& !< output file with scalar variables
        outputs_GR,& !< output file for the surrogate model fitting
        outputs_c,& !< output file with concentration profiles
        outputs_kin,& !< output file with variables of detailed kinetic model
        outputs_drain,& !< output file for the wall drainage
        outputs_af,& !< final foam properties for foam cond. and foam aging
        geometry !< geometry 3D=spherical, 2D=cylindrical
    logical :: &
        inertial_term,& !<include inertial term in equations (t/f)
        solcorr,& !<use solubility correction on bubble radius (t/f)
        surfactantPresent,& !<surfactant is present - pcsaft caluclation (t/f)
        gelpoint,& !<gel point reached (t/f)
        dilution,& !<use dilution effect for kinetics (t/f)
        shooting !<am I using shooting method (t/f)
    integer :: &
        fi1,fi2,fi3,fi4,fi5,& !<file indices
        integrator,& !<integrator. 1=dlsode,2=dlsodes
        int_meth,& !< stiff or no? See MF for ODEPACK
        p,& !<number of internal nodes
        maxts,& !<maximum inner time steps between t and t+h
        its,& !<number of outer integration outer time steps
        visc_model,& !<viscosity model. 1=constant,2=Castro and Macosko,3=modena
        rhop_model,& !<polymer density model. 1=constant,2=modena
        itens_model,& !<interfacial tension model. 1=constant,2=modena
        kin_model,& !<reaction kinetics model.
        ngas,& !<number of dissolved gases
        co2_pos,& !<carbon dioxide position
        fceq,& !<first concentration equation (index)
        fpeq,& !<first pressure equation (index)
        lpeq,& !<last pressure equation (index)
        req,& !<radius equation (index)
        teq,& !<temperature equation (index)
        xOHeq,& !<polyol conversion equation (index)
        xWeq !<water conversion equation (index)
    real(dp) :: &
        mshco,& !<mesh coarsening parameter
        temp0,& !<initial temperature
        R0,& !<initial radius
        nb0,& !<initial bubble number density
        Sn,& !<how many times is initial shell larger than initial bubble radius
        OH0,& !<initial concentration of polyol (don't set to zero)
        W0,& !<initial concentration of water (can be zero)
        NCO0,& !<initial concentration of isocyanate
        catalyst,& !<concentration of catalyst
        polyol1_ini,& !<initial concentration of polyol 1
        polyol2_ini,& !<initial concentration of polyol 2
        amine_ini,& !<initial concentration of amine
        isocyanate1_ini,& !<initial concentration of isocyanate 1
        isocyanate2_ini,& !<initial concentration of isocyanate 2
        isocyanate3_ini,& !<initial concentration of isocyanate 3
        AOH,& !<frequential factor of gelling reaction
        EOH,& !<activation energy of gelling reaction
        AW,& !<frequential factor of blowing reaction
        EW,& !<activation energy of blowing reaction
        dHOH,& !<gelling reaction enthalpy
        dHW,& !<blowing reaction enthalpy
        time,& !<time (s)
        radius,& !<bubble radius (m)
        laplace_pres,& !<Laplace pressure (Pa)
        eqconc,& !<equivalent concentration for first gas
        grrate(2),& !<growth rate
        st,& !<shell thickness
        S0,& !<initial shell thickness
        rel_tol,& !<relative tolerance
        abs_tol,& !<absolute tolerance
        eta,& !<viscosity
        pamb,& !<ambient pressure (in the liquid)
        sigma,& !<interfacial tension
        rhop,& !<polymer density
        cp,& !<heat capacity of the reaction mixture
        cppol,& !<heat capacity of polymer
        rhobl,& !<density of liquid physical blowing agent
        porosity,& !<foam porosity
        rhofoam,& !<foam density
        pair0,& !<initial partial pressure of air
        pair,& !<partial pressure of air
        timestep,& !<timestep
        goalRadius,& !<final radius we want to achieve
        nold(2),& !<moles in bubble
        vsh,& !<shell volume
        temp,& !<temperature (K)
        conv,& !<conversion of polyol
        gelpointconv,& !<conversion of polyol at gel point
        Rey,& !<Reynolds number
        pairst,& !<non-dimensional air pressure
        pambst,& !<non-dimensional ambient pressure
        Ca !<capillary number
    integer, dimension(:), allocatable :: &
        diff_model,& !<diffusivity model 1=constant,2=modena
        sol_model,& !<solubility model 1=constant,2=modena
        kineq !<kinetics state variable equations (indexes)
    real(dp), dimension(:), allocatable :: &
        y,& !<state variables
        cbl,& !<concentration profile in reaction mixture
        xgas,& !<initial molar fraction of gases in the bubble
        kinsource,& !<kinetic source term
        D,& !<diffusion coefficients (for each dissolved gas)
        D0,& !<initial diffusion coefficients used in non-dimensional routines
        KH,& !<Henry constants (for each dissolved gas)
        Mbl,& !<blowing agent molar mass (for each dissolved gas)
        dHv,& !<evaporation heat of blowing agent (for each gas)
        cpblg,& !<heat capacity of blowing agent in gas phase (for each gas)
        cpbll,& !<heat capacity of blowing agent in liquid phase (for each gas)
        mb,& !<moles in polymer
        mb2,& !<moles in bubble
        mb3,& !<total moles
        avconc,& !<average concentration in reaction mixture
        pressure,& !<partial pressure(Pa)
        wblpol,& !<weight fraction of blowing agents in reaction mixture
        dz !<spatial discretization
end module globals
