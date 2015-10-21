!>This file contains the subroutines to set the initial density profile.



 
!>The subroutine FormFunction is a wrapper function which takes care of 
!!handling the global PETSc data structures and creates local copys for 
!!every processor. It then calls the subroutine InitialGuessLocal which
!!sets the initial density values.
 
Subroutine FormInitialGuess(snes,X,ierr)
Use PetscManagement
Use f90moduleinterfaces

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
      type(userctx), pointer:: puser
      Vec            X
      PetscErrorCode ierr
      DM             da

!  Declarations for use with local arrays:
      PetscScalar,pointer :: lx_v(:,:)
      Vec                 :: localX

      ierr = 0
      call SNESGetDM(snes,da,ierr)
      call SNESGetApplicationContext(snes,puser,ierr)

!  Get a pointer to vector data.
      call DMGetLocalVector(da,localX,ierr)
      !call VecGetArrayF90(localX,lx_v,ierr) !only for DOF=1
      call DMDAVecGetArrayF90(da,localX,lx_v,ierr)


      
!  Compute initial guess over the locally owned part of the grid
      call InitialGuessLocal(puser,lx_v,ierr)

!  Restore vector
      !call VecRestoreArrayF90(localX,lx_v,ierr) !only for DOF=1
      call DMDAVecRestoreArrayF90(da,localX,lx_v,ierr)


!  Insert values into global vector
      call DMLocalToGlobalBegin(da,localX,INSERT_VALUES,X,ierr)
      call DMLocalToGlobalEnd(da,localX,INSERT_VALUES,X,ierr)
      call DMRestoreLocalVector(da,localX,ierr)

End Subroutine FormInitialGuess





!>In this subroutine every processor sets the initial density values
!!at its locally owned part of the grid.

Subroutine InitialGuessLocal(user,rhop,ierr)

!PETSc modules
      use PetscManagement
!VLE and DFT modules
      Use BASIC_VARIABLES, ONLY: ncomp,t,parame
      Use EOS_VARIABLES, Only:rho
      Use VLE_VAR, Only: rhob,tc
      Use PARAMETERS, Only: PI
      Use mod_DFT, Only: box,zp,dzp
      Use Global_x, Only: ngrid


      implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscdmda.h>
#include <finclude/petscis.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>
#include <finclude/petscsnes.h>

!  Input/output variables:
      type (userctx)   ::   user
      PetscScalar      ::   rhop(ncomp,user%gxs:user%gxe)
      PetscErrorCode   ::   ierr


!  Local variables:
      PetscInt  i,k
      PetscBool flg
      REAL    :: arg
      INTEGER :: pert
      REAL    :: tanhfac
      REAL    :: zp_i
      REAL    :: zp_middle



!normal or perturbed inital profile?
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-init_pert',pert,flg,ierr)


    zp_middle = (ngrid/2)*dzp        
    tanhfac   = -2.3625*t/tc + 2.4728


   !default: no perturbation, start with bulk density


IF(user%num_procs == 1) THEN !if only one processor involved
 
     Do k = 1,ncomp 
        Do i=user%gxs,user%gxe  !proc 0 has to calculate values for ghost points out of physical domain (from -irc to -1 and discret to discret+irc) and its regular part
              zp_i = zp(i)   
              rhop(k,i) =  ( TANH(-(zp_i -zp_middle) / parame(k,2) *tanhfac) + 1.0 ) * rhob(1,k)/2.0 &
                  - ( TANH(-(zp_i - zp_middle) / parame(k,2) *tanhfac) - 1.0 ) * rhob(2,k)/2.0
!rhop(k,i) = rhob(1,k)
              !make sure density is nonzero
              If(rhop(k,i) == 0.0) rhop(k,i) = rhop(k,i-1)  
        End Do
    End Do



Else !parallel simulation



    IF(user%rank == 0) THEN        !watch the user%xs vs user%gxs (first no ghost points, second with ghost points)      

        Do k=1,ncomp !loop over dof = loop over components
           Do i=user%gxs,user%xe  !proc 0 has to calculate values for ghost points out of physical domain (from -irc to -1) and its regular part
             zp_i = zp(i) 
             !perturb =  (rhob(1,j)-rhob(2,j))/2.0 * i * (discret - i) / (damp*discret)**2
             !rhop(k,i) =  ( TANH(-(zp_i -zp_middle) / parame(k,2) *tanhfac) + 1.0 ) * rhob(1,k)/2.0 &
             !    - ( TANH(-(zp_i - zp_middle) / parame(k,2) *tanhfac) - 1.0 ) * rhob(2,k)/2.0 !+ perturb
rhop(k,i) = rhob(1,k)
             !make sure density is nonzero
             If(rhop(k,i) == 0.0) rhop(k,i) = rhop(k,i-1)  
           End Do
        End Do 


    ELSE IF(user%rank == user%num_procs-1) THEN     

        Do k=1,ncomp !loop over dof = loop over components
           Do i=user%gxs,user%gxe  !last proc has to calculate values for ghost points out of physical domain (from discret to discret+irc) and its regular part
              zp_i = zp(i) 
              !perturb =  (rhob(1,j)-rhob(2,j))/2.0 * i * (discret - i) / (damp*discret)**2
              !rhop(k,i) =  ( TANH(-(zp_i -zp_middle) / parame(k,2) *tanhfac) + 1.0 ) * rhob(1,k)/2.0 &
              !   - ( TANH(-(zp_i - zp_middle) / parame(k,2) *tanhfac) - 1.0 ) * rhob(2,k)/2.0 !+ perturb
rhop(k,i) = rhob(1,k)
              !make sure density is nonzero
              If(rhop(k,i) == 0.0) rhop(k,i) = rhop(k,i-1)  
           End Do
        End Do 


     ELSE  

        Do k=1,ncomp !loop over dof = loop over components
           Do i=user%xs,user%xe  !middle procs have to calculate only their regular part
              zp_i = zp(i) 
              !perturb =  (rhob(1,j)-rhob(2,j))/2.0 * i * (discret - i) / (damp*discret)**2
              !rhop(k,i) =  ( TANH(-(zp_i -zp_middle) / parame(k,2) *tanhfac) + 1.0 ) * rhob(1,k)/2.0 &
              !   - ( TANH(-(zp_i - zp_middle) / parame(k,2) *tanhfac) - 1.0 ) * rhob(2,k)/2.0 !+ perturb
rhop(k,i) = rhob(1,k)
              !make sure density is nonzero
              If(rhop(k,i) == 0.0) rhop(k,i) = rhop(k,i-1)  
            End Do
        End Do 

     END IF


End If



















!       If(pert == 1) Then
!         Do i=user%xs,user%xe
!           arg = zp(i) * PI / box
!           rhop(k,i) = rhob(1,k) + 0.5*rhob(1,k) * sin(arg)
!         End Do
!       End If
!       
!       If(pert == 2) Then
!         Do i=user%xs,user%xe
!           arg = zp(i) * 2.0 * PI / box
!           rhop(k,i) = rhob(1,k) + 0.5*rhob(1,k) * sin(arg)
!         End Do
!       End If
!       
!       If(pert == 3) Then
!         Do i=user%xs,user%xe
!           arg = zp(i) * 3.0 * PI / box
!           rhop(k,i) = rhob(1,k) + 0.5*rhob(1,k) * sin(arg)
!         End Do
!       End If 




   







End Subroutine InitialGuessLocal



















