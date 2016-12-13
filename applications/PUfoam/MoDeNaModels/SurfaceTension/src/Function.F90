!> \file Function.F90
!! \brief This file contains the residual function of the density functional theory calculation.




!>The subroutine FormFunction is a wrapper function which takes care of 
!!handling the global PETSc data structures and creates local copys for 
!!every processor. It then calls the subroutine FormFunctionLocal which
!!performs the actual calculation of the residual at every grid point.

 
Subroutine FormFunction(snes,X,F,user,ierr)
Use PetscManagement
Implicit None
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>
#include <finclude/petscis.h>
#include <finclude/petscsnes.h>
#include <finclude/petscsnes.h90>

!  Input/output variables:
      SNES           snes
      Vec            X,F
      PetscErrorCode ierr
      type (userctx) user
      DM             da

!  Declarations for use with local arrays:
      PetscScalar,pointer :: lx_v(:,:),lf_v(:,:)
      Vec            localX

!  Scatter ghost points to local vector, using the 2-step process
!     DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
!  By placing code between these two statements, computations can
!  be done while messages are in transition.
      call SNESGetDM(snes,da,ierr)
      call DMGetLocalVector(da,localX,ierr)
      call DMGlobalToLocalBegin(da,X,INSERT_VALUES,                     &
     &     localX,ierr)
      call DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX,ierr)


      !call VecGetArrayF90(localX,lx_v,ierr)  !only for DOF=1
      !call VecGetArrayF90(F,lf_v,ierr)       !only for DOF=1
      call DMDAVecGetArrayF90(da,localX,lx_v,ierr)
      call DMDAVecGetArrayF90(da,F,lf_v,ierr)

!  Compute function over the locally owned part of the grid
      call FormFunctionLocal(lx_v,lf_v,user,ierr)

!  Restore vectors
      !call VecRestoreArrayF90(localX,lx_v,ierr) !only for DOF=1
      !call VecRestoreArrayF90(F,lf_v,ierr)      !only for DOF=1
      call DMDAVecRestoreArrayF90(da,localX,lx_v,ierr)
      call DMDAVecRestoreArrayF90(da,F,lf_v,ierr)

      
!  Insert values into global vector

      call DMRestoreLocalVector(da,localX,ierr)

      return
End Subroutine FormFunction











!>This subroutine performs the local evaluation of the residual function.
!!It calls the subroutines which calculate the different contributions to
!!the Helmholtz energy functional.


Subroutine FormFunctionLocal(rhop,f,user,ierr)
!PETSc modules
 Use PetscManagement
!DFT modules
 Use mod_DFT_FMT
 Use mod_DFT_CHAIN
 Use mod_DFT_DISP_WDA
 Use PARAMETERS, Only: PI ,muhs,muhc,mudisp
 Use BASIC_VARIABLES, Only: ncomp,nc,np  ,nphas,xi,dense,parame,ensemble_flag !nphas nur fur fugacity call 
 Use EOS_VARIABLES, Only:dhs,rho ,phas,eta_start,x,mseg,eta !letze 4 nur zum Test for aufruf eos_phi!!  
 Use mod_DFT, Only: fa,zp,dzp,free,pbulk, fa_disp, ab_disp !letzen beiden nur fur elmars version
 Use VLE_VAR, Only: rhob
 Use DFT_FCN_MODULE, Only: ChemPot_res,ChemPot_total
 Use Global_x, Only: ngrid, ngp

 
Implicit None

!  Input/output variables:
      type (userctx) user
      PetscScalar  rhop(ncomp,user%gxs:user%gxe)
      PetscScalar  f(ncomp,user%xs:user%xe)
      PetscErrorCode ierr

!  Local variables:
      PetscInt  i
      PetscInt  k

      !FMT
      REAL,dimension(user%gxs:user%gxe)          :: n0,n1,n2,n3,nv1,nv2  !ngp muss groesser als fa+fa/2 sein!!
      REAL,dimension(user%gxs:user%gxe)          :: phi_dn0,phi_dn1,phi_dn2,phi_dn3,phi_dnv1,phi_dnv2
      REAL                                       :: f_fmt,temp,zs  
      REAL,dimension(user%xs:user%xe,ncomp)      :: dF_drho_FMT

      !Chain
      REAL,dimension(user%gxs:user%gxe,ncomp)    :: rhobar,lambda  
      REAL                                       :: f_ch
      REAL,dimension(user%xs:user%xe,ncomp)      :: dF_drho_CHAIN

      !DISP VAR
      REAL,dimension(user%gxs:user%gxe,ncomp)    :: rhop_wd,my_disp,df_disp_drk
      REAL,dimension(user%gxs:user%gxe)          :: f_disp
      REAL, dimension(ncomp)                     :: dF_drho_disp
          !for Elmars Version
          REAL, DIMENSION(user%gxs:user%gxe,ncomp):: rhoi_disp,mydisp,dadisp_dr
          REAL, DIMENSION(user%gxs:user%gxe)     :: adisp,rho_disp,rhop_sum
          REAL                                   :: dF_drho_att(ncomp) 
          INTEGER                                :: fa_psi_max,WDA_var
          

      !polar 
      REAL                                       :: fres_polar,fdd,fqq,fdq,mu_polar(nc)
      REAL                                       :: fdd_rk, fqq_rk, fdq_rk, z3_rk  
      INTEGER                                    :: ik
      !association
      REAL                                       :: f_assoc
      REAL                                       :: mu_assoc(nc)
      REAL                                       :: lnphi(np,nc), fres
      
      REAL :: f_tot, delta_f      
      REAL :: Vext(ncomp)


      Do k = 1,ncomp
        Do i= user%xs,user%xe
              If( rhop(k,i) < epsilon(dhs) ) rhop(k,i) = epsilon(dhs)
        End Do
      End Do

      DO i = user%gxs,user%gxe
      rhop_sum(i) = sum( rhop(1:ncomp,i) )
     END DO



      !calculate weighted densities
      call FMT_Weighted_Densities(rhop,n0,n1,n2,n3,nv1,nv2,phi_dn0,phi_dn1,phi_dn2,phi_dn3,phi_dnv1,phi_dnv2,user)
      
      !calculate averaged density rhobar and lambda (both needed for chain term)
      call Chain_aux(rhop,rhobar,lambda,user) 

        
      !weighted densities for Dispersion 
!       call DISP_Weighted_Densities(rhop, rhop_wd, user)
!       call DISP_mu(rhop_wd,f_disp,my_disp,df_disp_drk,user)
        
        !Elmars Version 
        fa_psi_max = maxval(fa_disp(1:ncomp))   
        call rhoi_disp_wd ( ngrid, fa_disp, fa_psi_max, ab_disp, rhop, rhoi_disp,user )
        CALL a_disp_pcsaft ( ngrid, fa_disp, fa_psi_max, rhoi_disp, rho_disp, adisp, mydisp, dadisp_dr,user )

      
      free = 0.0
      
      DO i = user%xs,user%xe

           !FMT
           call FMT_dFdrho(i,fa,user,phi_dn0,phi_dn1,phi_dn2,phi_dn3,phi_dnv1,phi_dnv2,dF_drho_FMT)
           zs = 1.0 - n3(i)
           temp = log(zs)           
           f_fmt = - n0(i)*temp + n1(i)*n2(i)/zs - nv1(i)*nv2(i)/zs  &              !see for example eq. 30 in "Fundamental measure theory for hard-sphere mixtures revisited: the White bear version (Roth)
            + (n2(i)**3 -3.0*n2(i)*nv2(i)*nv2(i)) *(n3(i)+zs*zs*temp)  &
            /36.0/PI/zs/zs/n3(i)**2 

           !Chain
           call Chain_dFdrho(i,rhop,lambda,rhobar,dF_drho_CHAIN,f_ch,user)

           !Dispersion
           WDA_var = 1
           call dF_disp_drho_wda( i, WDA_var, fa_disp, ab_disp, rhop, rhop_sum, &
                         rhoi_disp, rho_disp, adisp, mydisp, dadisp_dr, dF_drho_att,user )           

           if(WDA_var == 1) adisp(i) = rho_disp(i) * adisp(i)   
           if(WDA_var == 2) adisp(i) = rhop_sum(i) * adisp(i)   


           !polar as LDA
           fres_polar  = 0.0
           mu_polar(:) = 0.0 
           dense(1) = PI / 6.0 * SUM( rhop(1:ncomp,i) * mseg(1:ncomp) * dhs(1:ncomp)**3 )
           xi(1,1:ncomp) = rhop(1:ncomp,i) / SUM( rhop(1:ncomp,i) )    
           ensemble_flag = 'tv' 
           IF ( SUM( parame(1:ncomp,6) ) > 1.E-10 .OR. SUM( parame(1:ncomp,7) ) > 1.E-10 ) THEN
              eta = dense(1) 
              rho = SUM(rhop(1:ncomp,i))
              x(1:ncomp) = xi(1,1:ncomp)
              call F_POLAR ( fdd, fqq, fdq )
              fres_polar  = ( fdd + fqq + fdq ) * SUM( rhop(1:ncomp,i) )
              DO ik = 1, ncomp
                  z3_rk = PI/6.0 * mseg(ik) * dhs(ik)**3 
                  call PHI_POLAR ( ik, z3_rk, fdd_rk, fqq_rk, fdq_rk )
                   mu_polar(ik) = fdd_rk + fqq_rk + fdq_rk
              END DO
            END IF

           !association as LDA
           f_assoc  = 0.0
           mu_assoc(:) = 0.0
           IF ( SUM( NINT(parame(1:ncomp,12)) ) > 0) THEN
              call ONLY_ONE_TERM_EOS_NUMERICAL ( 'hb_term  ', 'TPT1_Chap' )
              call FUGACITY ( lnphi )
              call RESTORE_PREVIOUS_EOS_NUMERICAL
              f_assoc  = fres * SUM( rhop(1:ncomp,i) )
              mu_assoc(1:ncomp) = lnphi(1,1:ncomp)
           END IF      

           f_tot = f_fmt + f_ch + adisp(i) + fres_polar + f_assoc &
                         + SUM( rhop(1:ncomp,i)*( LOG(rhop(1:ncomp,i)/rhob(1,1:ncomp))-1.0 ) )
           delta_f = f_tot   -   ( SUM(rhop(1:ncomp,i)*chemPot_res(1:ncomp)) - pbulk)  ! all quantities .../(kT)
           free  = free  + delta_f*dzp           



           Do k=1,ncomp
               f(k,i) = -(rhob(1,k) * exp(ChemPot_res(k)-dF_drho_FMT(i,k)-dF_drho_CHAIN(i,k)-dF_drho_att(k) &
                                            -mu_polar(k) - mu_assoc(k)) - rhop(k,i))
           End Do

       END DO



End Subroutine FormFunctionLocal






