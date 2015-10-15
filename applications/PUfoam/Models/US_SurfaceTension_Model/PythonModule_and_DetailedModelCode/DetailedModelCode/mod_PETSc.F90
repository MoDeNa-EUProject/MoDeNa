



!>In this module, the application context is defined

Module PetscManagement

      type userctx
#include <finclude/petscsysdef.h>
#include <finclude/petscvecdef.h>
#include <finclude/petscdmdef.h>
        PetscInt xs,xe,xm,gxs,gxe,gxm,mx
        PetscMPIInt rank
        PetscMPIInt num_procs
      end type userctx

End Module
 


!>This module contains variables associated with the PETSc solver

Module Global_x
Implicit None
#include <finclude/petscsysdef.h>
#include <finclude/petscsnesdef.h>
#include <finclude/petscvecdef.h>    

SNES     :: snes                 !nonlinear solver context
Vec      :: x                    !array of unknowns  
Vec      :: r                    !array of residuals
PetscInt :: ngrid                !number of grid points
PetscInt :: ngp                  !number of ghost points
 
DOUBLE PRECISION :: timer        ! timer
DOUBLE PRECISION :: timer_old    ! timer
DOUBLE PRECISION :: total_time   ! timer


End Module Global_x



Module f90moduleinterfaces
Use PetscManagement
#include <finclude/petscsnesdef.h>

      Interface SNESSetApplicationContext
        Subroutine SNESSetApplicationContext(snes,ctx,ierr)
        use PetscManagement
          SNES snes
          type(userctx) ctx
          PetscErrorCode ierr
        End Subroutine
      End Interface SNESSetApplicationContext

      Interface SNESGetApplicationContext
        Subroutine SNESGetApplicationContext(snes,ctx,ierr)
        use PetscManagement
          SNES snes
          type(userctx), pointer :: ctx
          PetscErrorCode ierr
        End Subroutine
      End Interface SNESGetApplicationContext

End Module f90moduleinterfaces


