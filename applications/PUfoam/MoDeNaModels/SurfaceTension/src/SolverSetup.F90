

!> \file SolverSetup.f90 
!!This file contains the subroutine which sets the method used to
!!solve the nonlinear system of equations depending on the specifications
!!made in the makefile.



 
Subroutine SolverSetup

Use BASIC_VARIABLES, Only: ncomp

!PETSc modules
Use PetscManagement
Use f90moduleinterfaces
Use Global_x      !(snes,ngrid,ngp,x,r,timer)

Implicit None

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscdmda.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h90>
#include <finclude/petscis.h>
#include <finclude/petscmat.h>
#include <finclude/petscsnes.h>
#include <finclude/petscsnes.h90>


DM             ::  da
PetscErrorCode ::  ierr
type (userctx) ::  user
Mat            ::  J
PetscReal      ::  erel
PetscBool      ::  flg
INTEGER        ::  jac_id
PetscInt       ::  JacLocal

character(80) :: filename=''

!external subroutines associated with Solver
external FormInitialGuess
external FormFunction
external Jac_Shell_AD
external Jac_Matrix_Empty
external MonitorTimer


     call MPI_Comm_rank(PETSC_COMM_WORLD,user%rank,ierr)
     call MPI_Comm_size(PETSC_COMM_WORLD,user%num_procs,ierr) 


!  Create solver context
     call SNESCreate(PETSC_COMM_WORLD,snes,ierr) 

!  Create distributed array (DMDA) to manage parallel grid and vectors
     call DMDACreate1d(PETSC_COMM_WORLD,                                &  !MPI communicator
     &     DMDA_BOUNDARY_GHOSTED,                                          &  !Boundary type at boundary of physical domain
     &     ngrid,                                                       &  !global dimension of array (if negative number, then value can be changed by user via command line!)
     &     ncomp,                                                           &  !number of degrees of freedom per grid point (number of unknowns at each grid point)
     &     ngp,                                                           &  !number of ghost points accessible for local vectors
     &     PETSC_NULL_INTEGER,                                          &  !could be an array to specify the number of grid points per processor
     &     da,                                                          &  !the resulting distributed array object 
     &     ierr)   

!  Extract global arrays from DMDA: x: unknowns, r: residual
      call DMCreateGlobalVector(da,x,ierr)
      call VecDuplicate(x,r,ierr)

!  Get local grid boundaries (for 1-dimensional DMDA)
 call DMDAGetCorners(da,       & !the distributed array 
     &   user%xs,                   & !corner index in x direction
     &   PETSC_NULL_INTEGER,        & !corner index in y direction
     &   PETSC_NULL_INTEGER,        & !corner index in z direction
     &   user%xm,                   & !width of locally owned part in x direction
     &   PETSC_NULL_INTEGER,        & !width of locally owned part in y direction
     &   PETSC_NULL_INTEGER,        & !width of locally owned part in z direction 
     &   ierr)                        !error check 
       
      call DMDAGetGhostCorners(da,  & !the distributed array
     &   user%gxs,                  & !corner index in x direction (but now counting includes ghost points) 
     &   PETSC_NULL_INTEGER,        & !corner index in y direction (but now counting includes ghost points)
     &   PETSC_NULL_INTEGER,        & !corner index in z direction (but now counting includes ghost points)
     &   user%gxm,                  & !width of locally owned part in x direction (but now including ghost points)
     &   PETSC_NULL_INTEGER,        & !width of locally owned part in y direction (but now including ghost points)
     &   PETSC_NULL_INTEGER,        & !width of locally owned part in z direction (but now including ghost points)
     &   ierr)                        !error check

!  Here we shift the starting indices up by one so that we can easily
!  use the Fortran convention of 1-based indices (rather 0-based indices).
      user%xs  = user%xs+1
      user%gxs = user%gxs+1
      user%xe  = user%xs+user%xm-1
      user%gxe = user%gxs+user%gxm-1

      call SNESSetFunction(snes,r,FormFunction,user,ierr)
      call SNESSetApplicationContext(snes,user,ierr)
      call SNESSetDM(snes,da,ierr)

! !Set up matrix free jacobian
!       erel = 1e-08
!       call MatCreateSNESMF(snes,J,ierr)  !matrix free jacobi matrix
!       call MatMFFDSetFunctionError(J,erel,ierr)
!       call SNESSetJacobian(snes,J,J,MatMFFDComputeJacobian,PETSC_NULL_OBJECT,ierr)


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Set up nonlinear solver depending on option set in makefile (with -jac)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-jac',jac_id,flg,ierr)
      If(flg) Then
                
      Else         
             jac_id = 0
             write(*,*)'Option -jac not set in makefile, using matrix-free with numerical approximations.'
      End If
      
      Select Case(jac_id)
       
            Case(-1) !for options that dont need a jacobian: Anderson mixing, Picard

            Case (0) !matrix-free with numerical approximation of J(x)v        
               !default value of erel
               erel = 1e-08
               !check if value for erel is specified in makefile
               call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-erel',erel,flg,ierr)      
               If(user%rank == 0) write(*,*)'Using matrix-free Jacobian with numerical approximations. erel:',erel
               call MatCreateSNESMF(snes,J,ierr)  !matrix free jacobi matrix
               call MatMFFDSetFunctionError(J,erel,ierr)
               call SNESSetJacobian(snes,J,J,MatMFFDComputeJacobian,PETSC_NULL_OBJECT,ierr)
 
            Case (1) !matrix-free with AD-calculated J(x)v
               If(user%rank == 0) write(*,*)'Using matrix-free AD Jacobian.'
               !determine local size of shell matrix
               JacLocal = INT(ngrid / user%num_procs) 
               If(mod(ngrid,user%num_procs) /= 0 ) Then
                   write(*,*)'Surface Tension Code: Dimensions and number of cores dont match! mod(ngrid,n_cores) must equal 0'
                   Stop 5                        
               End if
               
               !Make sure that local parts of JacShell have the same size as the local DMDA-Arrays (local size JacShell != user%xm!!)
               If(JacLocal /= user%xm) Then
                   write(*,*)'Surface Tension Code: Shell-Jacobi-Matrix and DMDA have to be parallelized accordingly.'
!               write(*,*) 'Surface Tension Code: Shell-Jacobi-Matrix and DMDA have to be parallelized accordingly!!(JacLocal = user%xm)'
                   Stop 6 
               End If
               
               call MatCreateShell(PETSC_COMM_WORLD,ncomp*JacLocal,ncomp*JacLocal,ncomp*ngrid,ncomp*ngrid,PETSC_NULL_OBJECT,J,ierr)
               call MatShellSetOperation( J, MATOP_MULT, Jac_Shell_AD, ierr )
               call SNESSetJacobian( snes, J, J, Jac_Matrix_Empty, PETSC_NULL_OBJECT, ierr)
              
           Case (2) !build complete Jacobi matrix via finite-differences
              If(user%rank == 0) write(*,*)'Using finite-difference Jacobian.'
!Now in makefile: -snes_fd
!               call MatCreate(PETSC_COMM_WORLD,J,ierr)
!               call MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,ngrid,ngrid,ierr)  
!               call MatSetFromOptions(J,ierr)                            
!               call MatSetUp(J,ierr)      !!sets up the internal matrix data structure 
!               call SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault,PETSC_NULL_OBJECT,ierr)

               
           Case (3) !build the complete Jacobi matrix via AD
               If(user%rank == 0) write(*,*)'Surface Tension Code: Using AD-generated Jacobian.'
               write(*,*)'Not working yet!' 
               stop 6              
               !create Matrix that has the appropriate sparsity pattern for the da
               !But that also means that when values are inserted in positions that
               !dont fit this structure, an error occurs!
               !i.e. when calculating the jacobi matrix it is not possible anymore
               !to ignore the sparsity pattern and just create a dense jacobian!!                 
               
               !call DMCreateMatrix(da,MATMPIAIJ,J,ierr)   !CAUTION: in newer PETSc Versions DMCreateMatrix does not contain the MatType argument anymore!
               !call SNESSetJacobian(snes,J,J,Jac_AD,PETSC_NULL_OBJECT,ierr)
            
          Case default
              write(*,*) 'No valid option set for -jac in makefile, using full Jacobi via finite-differences.'  
              call SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault,PETSC_NULL_OBJECT,ierr)
 
      End Select

!Set monitoring function that ouputs the surface tension at every iteration
      call SNESMonitorSet(snes,MonitorTimer,PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)

!Enable modification from makefile
      call SNESSetFromOptions(snes,ierr)

!Evaluate Initial Guess
      call FormInitialGuess(snes,x,ierr)
      write(filename,*)'0_initial_profile_global.xlo'
      call PrintGlobalVec(x,filename)

      !start timer
      total_time = 0.0
      timer_old = MPI_WTIME()

!Solve system
      call SNESSolve(snes,PETSC_NULL_OBJECT,x,ierr)
      write(filename,*)'1_final_profile_global.xlo'
      call PrintGlobalVec(x,filename)
      write(filename,'(a,I3.3,a)') '2_final_profile_local_proc_',user%rank,'.xlo' 
      call PrintLocalVec(x,da,filename)  

      !plot residual vs time graph
      !call system('gnuplot gnuplot_script.srp')  
    

End Subroutine SolverSetup







subroutine MonitorTimer(snes,its,norm,dummy,ierr)

  Use Global_x, Only: timer,timer_old,total_time,r
  Use mod_DFT,  Only: free    
  Use PARAMETERS, Only: KBOL,PI
  Use BASIC_VARIABLES, Only: t,parame, ncomp
  Use VLE_VAR, Only: rhob,tc 

! ! !PETSc modules
! !       use f90module   
! ! !DFT modules
! !       use DFT_MODULE, ONLY: free
! !       use VLE_VAR, ONLY: tc,rhob
! !       use BASIC_VARIABLES, ONLY:t,parame,ncomp
! !       use PARAMETERS, ONLY: PI,KBOL   

      implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>

!  Input/output variables:
      SNES           snes
      PetscInt   dummy
      Integer :: its
      REAL    :: norm
      PetscErrorCode ierr
! local
      PetscReal :: f2norm
      Vec       :: current_solution
      DOUBLE PRECISION :: delta_time
      INTEGER   :: rank
      REAL      :: m_average,surftens,st_macro
      character(80) :: filename=''



        !calculate interfacial tension 
          !sum up the variable free from all procs on proc 0
          call MPI_Reduce(free,free,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,PETSC_COMM_WORLD,ierr)  

          !Only proc 0 calculates surface tension
          call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
          IF(rank == 0) THEN
    
             !calculate surface tension  
             surftens = KBOL * t *1.E20*1000.0 *free 
             m_average = SUM( rhob(1,1:ncomp)*parame(1:ncomp,1) ) / rhob(1,0)

             st_macro = surftens / ( 1.0 + 3.0/8.0/PI *t/tc  &
                                   * (1.0/2.55)**2  / (0.0674*m_average+0.0045) )
             write(*,*)'ST',st_macro

            !write result to outputfile
            filename='./out.txt'
            CALL file_open(filename,99)
            WRITE (99,*) st_macro
            close(99)
             
          End If
        

        !calculate elapsed time
        timer = MPI_WTIME()  
        delta_time = timer - timer_old
        timer_old = timer
        total_time = total_time + delta_time
        
        !get norm of residual array
        call VecNorm(r,2,f2norm,ierr)

        IF(rank == 0) THEN
           !print results to file
           filename = 'ItsTimeNorm.dat'
           open(unit = 44, file = filename)
           write(44,*) its,total_time,f2norm
        End If 
           
end subroutine MonitorTimer

















