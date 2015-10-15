
!>This file contains subroutines that make printing of 
!!local and global PETSc arrays easier.


 
Subroutine PrintGlobalVec(GlobalVec,filename)
Implicit None
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscis.h>


     Vec           :: GlobalVec
     character(80) :: filename

     PetscViewer  viewer
     PetscErrorCode ierr

!Print global vector to a file
     call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,viewer,ierr)
     call VecView(GlobalVec,viewer,ierr)

     call PetscViewerDestroy(viewer,ierr)

End Subroutine PrintGlobalVec



Subroutine PrintLocalVec(GlobalVec,da,filename)
Implicit None
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>
#include <finclude/petscis.h>

        Vec GlobalVec
        Vec LocalVec
        DM da
        character(80) :: filename

        PetscViewer viewer
        PetscErrorCode ierr

        call DMGetLocalVector(da,LocalVec,ierr)                             !create a locally owned part of a global distribued array
        call DMGlobalToLocalBegin(da,GlobalVec,INSERT_VALUES,LocalVec,ierr) !copy values from global array to local arrays
        call DMGlobalToLocalEnd(da,GlobalVec,INSERT_VALUES,LocalVec,ierr) 
        call PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,viewer,ierr)       !associate viewer with file
        call VecView(LocalVec,viewer,ierr)                                    !every processor prints his local array to a file
        call DMRestoreLocalVector(da,LocalVec,ierr)                           !Free memory of local vector
 

End Subroutine PrintLocalVec
