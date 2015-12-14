!> This file contains the subroutine which starts the phase equilibrium calculation.
!! It also prints the calculated density to the outputfile "out.txt".



SUBROUTINE VLE_MIX(rhob,density,chemPot_total,user)


!Petsc modules
 USE PetscManagement


!VLE modules
 USE parameters, ONLY: PI, RGAS, KBOL, muhs,muhc,mudisp
 USE basic_variables
 USE EOS_VARIABLES, ONLY: fres, eta, eta_start, dhs, mseg, uij, sig_ij, rho, x, z3t
 USE DFT_MODULE
 USE EOS_NUMERICAL_DERIVATIVES
 USE DFT_FCN_MODULE, ONLY: chemPot_res
 IMPLICIT NONE

!> ---------------------------------------------------------------------
!! Variables
!! ---------------------------------------------------------------------
  
!passed
 type (userctx)         user 
 REAL                                   :: chemPot_total(nc)
 REAL                                   :: rhob(2,0:nc),density(np)

!local
 REAL, DIMENSION(nc)                    :: dhs_star
 REAL                                   :: w(np,nc), lnphi(np,nc)
 INTEGER                                :: converg

 !> ---------------------------------------------------------------------
 !! prepare for phase equilibrium calculation for given T
 !! ---------------------------------------------------------------------

  dhs(1:ncomp) = parame(1:ncomp,2) * ( 1.0 - 0.12*EXP( -3.0*parame(1:ncomp,3)/t ) )  ! needed for rdf_matrix
  dhs_star(1:ncomp) = dhs(1:ncomp)/parame(1:ncomp,2)

  nphas  = 2
  outp = 0                      ! output to terminal

  CALL START_VAR (converg,user)      ! gets starting values, sets "val_init"

  IF ( converg /= 1 ) THEN
    IF(user%rank == 0) THEN
    WRITE (*,*) 'no VLE found'
    END IF
    RETURN
  END IF

! rhob(phase,0): molecular density
 rhob(1,0) = dense(1) / (  PI/6.0* SUM( xi(1,1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )  )
 rhob(2,0) = dense(2) / (  PI/6.0* SUM( xi(2,1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )  )
 ! rhob(phase,i): molecular component density (with i=(1,...ncomp) ) in units (1/A^3)
 rhob(1,1:ncomp) = rhob(1,0)*xi(1,1:ncomp)
 rhob(2,1:ncomp) = rhob(2,0)*xi(2,1:ncomp)

! --- get density in SI-units (kg/m**3) -------------------------------
 CALL SI_DENS ( density, w )

!--- calculate residual chemical potentials
 ensemble_flag = 'tv'                    ! this flag is for: mu_res=mu_res(T,rho)
 densta(1) = dense(1)                    ! Index 1 is for liquid density (here: packing fraction eta)
 densta(2) = dense(2)                    ! Index 2 is for vapour density (here: packing fraction eta)
 CALL fugacity (lnphi)
 chemPot_res(1:ncomp) = lnphi(1,1:ncomp)
 chemPot_total(1:ncomp) = lnphi(1,1:ncomp) + LOG( rhob(1,1:ncomp) )     ! my0 = mu_res(T,rho_bulk_L) + ln(rho_bulk_l)


!> ---------------------------------------------------------------------
!! Output results of phase equilibrim calculation
!! ---------------------------------------------------------------------
 
 
IF(user%rank == 0) THEN
 WRITE(*,*) '--------------------------------------------------'
 WRITE(*,*)'RESULT OF PHASE EQUILIBRIUM CALCULATION'
 WRITE (*,*) ' '
 WRITE (*,*) 'temperature  ',t,      'K,  and  p=', p/1.E5,' bar'
 WRITE (*,*) 'x1_liquid    ',xi(1,1),'   x1_vapor', xi(2,1)
 WRITE (*,*) 'densities    ',rhob(1,0), rhob(2,0)
 WRITE (*,*) 'dense        ',dense(1), dense(2)
 WRITE (*,*) 'density [kg/m3]        ',density(1), density(2)
 write (*,*) 'chemical potentials comp1' , lnphi(1,1) + LOG( rhob(1,1) ),  lnphi(2,1) +  lnx(2,1) + LOG(rhob(2,0)) !LOG( rhob(2,1) )  
 write (*,*) 'chemical potentials comp2' ,lnphi(1,2) + LOG( rhob(1,2) ), lnphi(2,2) + LOG( rhob(2,2) )     
END IF


                             

END SUBROUTINE VLE_MIX 