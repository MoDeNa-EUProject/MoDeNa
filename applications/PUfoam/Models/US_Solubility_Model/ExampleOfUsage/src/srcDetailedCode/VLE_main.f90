SUBROUTINE VLE_MIX(rhob,density,chemPot_total,compID)

 USE parameters, ONLY: PI, RGAS, KBOL
 USE basic_variables
 USE EOS_VARIABLES, ONLY: fres, eta, eta_start, dhs, mseg, uij, sig_ij, rho, x, z3t
 USE DFT_MODULE
 USE EOS_NUMERICAL_DERIVATIVES
 IMPLICIT NONE
 


! ---------------------------------------------------------------------
! Variables
! ---------------------------------------------------------------------
  
 !passed
 REAL                                   :: chemPot_total(nc)
 REAL                                   :: rhob(2,0:nc),density(np)
 INTEGER                                :: compID
 
 !local
 REAL, DIMENSION(nc)                    :: dhs_star
 REAL                                   :: w(np,nc), lnphi(np,nc)
 INTEGER                                :: converg
 CHARACTER(LEN=4)                       :: char_ncomp 
 REAL                                   :: Henry
 INTEGER                                :: i
 CHARACTER (LEN=50)                     :: filename
 

! ---------------------------------------------------------------------
! prepare for phase equilibrium calculation for given T
! ---------------------------------------------------------------------

    dhs(1:ncomp) = parame(1:ncomp,2) * ( 1.0 - 0.12*EXP( -3.0*parame(1:ncomp,3)/t ) )  ! needed for rdf_matrix
    dhs_star(1:ncomp) = dhs(1:ncomp)/parame(1:ncomp,2)

    nphas  = 2                      ! number of phases
    outp   = 0                      ! output to terminal

    CALL START_VAR (converg)         

    IF ( converg /= 1 ) THEN
      WRITE (*,*) 'no VLE found'
      RETURN
    END IF

   ! rhob(phase,0): molecular density
    rhob(1,0) = dense(1) / (  PI/6.0* SUM( xi(1,1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )  )
    rhob(2,0) = dense(2) / (  PI/6.0* SUM( xi(2,1:ncomp) * parame(1:ncomp,1) * dhs(1:ncomp)**3 )  )
   ! rhob(phase,i): molecular component density (with i=(1,...ncomp) ) in units (1/A^3)
    rhob(1,1:ncomp) = rhob(1,0)*xi(1,1:ncomp)
    rhob(2,1:ncomp) = rhob(2,0)*xi(2,1:ncomp)

   ! get density in SI-units (kg/m**3) 
    CALL SI_DENS ( density, w )

   ! calculate residual chemical potentials
    ensemble_flag = 'tv'                    ! this flag is for: mu_res=mu_res(T,rho)
    densta(1) = dense(1)                    ! Index 1 is for liquid density (here: packing fraction eta)
    densta(2) = dense(2)                    ! Index 2 is for vapour density (here: packing fraction eta)
    CALL fugacity (lnphi)                   ! calculate fugacities
 
 
! ---------------------------------------------------------------------
! Output results of phase equilibrim calculation
! ---------------------------------------------------------------------
 
   WRITE(*,*) '--------------------------------------------------'
   WRITE(*,*)'RESULT OF PHASE EQUILIBRIUM CALCULATION'
   WRITE (char_ncomp,'(I3)') ncomp
   WRITE (*,*) 'T = ',t,      'K,  and  p =', p/1.E5,' bar'
   WRITE(*,*)' '
   WRITE(*,'(t15,4(a12,1x),10x,a)') (compna(i),i=1,ncomp)
   !!WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE I  w', (w(1,i),i=1,ncomp)
   !!WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE II w', (w(2,i),i=1,ncomp)
   WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE I  x', (EXP(lnx(1,i)),i=1,ncomp)
   WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE II x', (EXP(lnx(2,i)),i=1,ncomp)
   WRITE(*,*)' '
   !!WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE I  density', (rhob(1,i),i=1,ncomp)
   !! WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE II density', (rhob(2,i),i=1,ncomp)
   !!WRITE(*,*)' '
   WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE I  chemPot', (lnphi(1,i) + LOG(rhob(1,i)),i=1,ncomp)
   WRITE(*,'(2x,a,'//char_ncomp//'(g13.6,1x))') 'PHASE II chemPot', (lnphi(2,i) + LOG(rhob(2,i)),i=1,ncomp)
   WRITE(*,*)' '
   !!WRITE(*,*)'Phase densities '
   WRITE(*,'(2x,a,2(g13.6,1x))')                'SI-DENSITY [kg/m3]   ', density(1),density(2)
   !!WRITE(*,'(2x,a,2(g13.6,1x))')                'NUMBER-DENSITY       ', rhob(1,0),rhob(2,0) 
   WRITE(*,*)

  !Calculate Solubility
  !output: Henrys law: xi*H(T) = yi * p -> H = yip/xi
   Henry = EXP(lnx(2,1))*p /EXP(lnx(1,1)) 
   
   !write solubility to outputfile "out.txt"
   filename='./out.txt'
   CALL file_open(filename,78)
   write(78,*) Henry / 100000.0
   write(*,*)'Henry coefficient (bar):',Henry / 100000.0
   ! write(*,*)'g_CO2 / g_PU', exp(lnx(1,2)) * mm(2) / (exp(lnx(1,1)) * mm(1))
  
END SUBROUTINE VLE_MIX 
