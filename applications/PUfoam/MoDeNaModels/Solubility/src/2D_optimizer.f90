MODULE optimizer_2D

  USE BASIC_VARIABLES, only: outp
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Newton_Opt_2D, Newton_Opt_analytic, f_grad_hessian, f_grad

CONTAINS


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine Newton_Opt_analytic ( f_grad_hessian, n, x, min_grad, eps, g, f, convergence )

    implicit none

    !-----------------------------------------------------------------------------
    integer, intent(IN)                    :: n
    real, intent(IN OUT)                   :: x(:)
    real, intent(IN OUT)                   :: f
    real, intent(IN OUT)                   :: g(:)
    real, intent(IN)                       :: min_grad
    real, intent(IN)                       :: eps
    integer, intent(OUT)                   :: convergence

    !---------------------------------------------------------------------------
    INTERFACE
       SUBROUTINE f_grad_hessian (n, x, f, g, H, diagonal)
         INTEGER, INTENT(IN)                 :: n
         REAL, INTENT(IN)                    :: x(n)
         REAL, INTENT(OUT)                   :: f
         REAL, INTENT(OUT)                   :: g(n)
         REAL, INTENT(OUT)                   :: H(n,n)
         REAL, INTENT(OUT)                   :: diagonal(n,n)
       END SUBROUTINE f_grad_hessian
    END INTERFACE

    !---------------------------------------------------------------------------
    integer                                 :: count, max_c
    integer                                 :: it_num, rot_num
    real                                    :: determinant
    real                                    :: error_in_x, error_in_g
    real, dimension(n)                      :: x0, delta_x, eigenval
    real, dimension(n,n)                    :: H, H_x, diagonal, eigenvectors
    logical                                 :: converged
    !-----------------------------------------------------------------------------

    convergence = 0
    max_c = 50                    ! maximum number of iterations
    count = 0
    converged = .false.

    do while ( .NOT.converged .AND. count <= max_c )

       count = count + 1

       x0(:) = x(:)

       call f_grad_hessian (n, x, f, g, H, diagonal)

       error_in_g = maxval( abs( g(1:n) ) )

       call jacobi_eigenvalue ( n, H, 100, eigenvectors, eigenval, it_num, rot_num )

       if ( eigenval(1) > 1.E-3 .OR. (eigenval(1) > 0.0 .AND. error_in_g < 1.E-2) ) then

          delta_x(:) = g(:)
          call MATINV( n, 1, H, delta_x, determinant )  ! output delta_x := inv( H ) * delta_x
          ! write (*,*) 'determinant', determinant
          ! write (*,*) 'delta_x', delta_x(:)

          x(:) = x0(:) - delta_x(:)

       else

          H_x = H + 2.0 * eigenval(1) * diagonal
          call jacobi_eigenvalue ( n, H_x, 100, eigenvectors, eigenval, it_num, rot_num )

          delta_x(:) = g(:)
          call MATINV( n, 1, H_x, delta_x, determinant )  ! output delta_x := inv( H ) * delta_x

          x(:) = x0(:) - delta_x(:)

       end if

       error_in_x = sum( abs( delta_x(1:n) ) )

       if ( error_in_x < eps .OR. error_in_g < min_grad ) then
          converged = .true.
          convergence = 1
       end if

       if (outp >= 1) write (*,'(a,4(G15.7))') ' finished Newton step: f, errors ',f, error_in_x, error_in_g

    end do


  end subroutine Newton_Opt_analytic



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  subroutine Newton_Opt_2D ( f_val, x, n, min_grad, eps, g, f )

    implicit none

    !-----------------------------------------------------------------------------
    integer, intent(IN)                    :: n
    real, intent(IN OUT)                   :: x(:)
    real, intent(IN OUT)                   :: f
    real, intent(IN OUT)                   :: g(:)
    real, intent(IN)                       :: min_grad
    real, intent(IN)                       :: eps

    !---------------------------------------------------------------------------
    INTERFACE
       !SUBROUTINE f_grad_hessian (H, g, f, x, n)
       !  REAL, INTENT(IN OUT)                :: H(:,:)
       !  REAL, INTENT(IN OUT)                :: g(:)
       !  REAL, INTENT(IN OUT)                :: f
       !  REAL, INTENT(IN)                    :: x(:)
       !  INTEGER, INTENT(IN)                 :: n
       !END SUBROUTINE f_grad_hessian
       SUBROUTINE f_val (f, x, n)
         REAL, INTENT(IN OUT)                :: f
         REAL, INTENT(IN)                    :: x(:)
         INTEGER, INTENT(IN)                 :: n
       END SUBROUTINE f_val
    END INTERFACE

    !---------------------------------------------------------------------------
    integer                                 :: i, count, max_c
    integer                                 :: line_steps
    real                                    :: ftrial, det
    real                                    :: error_in_x
    real, dimension(n)                      :: x0
    real, dimension(n,n)                    :: H, invers_H
    real                                    :: ft, ft_min, length_scale
    real, dimension(n)                      :: xt, direction
    logical                                 :: line_search
    !-----------------------------------------------------------------------------

    if (n /= 2) write (*,*) 'Newton_Opt_2D is only defined for n=2'
    if (n /= 2) stop

    max_c = 50                    ! maximum number of iterations
    count = 0
    error_in_x = 10.0 * eps
    line_steps = 5

    do while ( error_in_x > eps .AND. count <= max_c )

       count = count + 1

       x0 = x

       call f_grad_hessian ( H, f_val, g, f, x, n )

       det = H(1,1) * H(2,2) - H(1,2) * H(2,1)
       if ( SUM( ABS( g(1:n) ) ) < min_grad .OR. ABS(det) < 1.E-20) EXIT

       invers_H(1,1) =   H(2,2)  / det
       invers_H(2,2) =   H(1,1)  / det
       invers_H(1,2) = - H(2,1)  / det
       invers_H(2,1) = - H(1,2)  / det

       if ( det > 0.0 ) then
          do i = 1,n
             x(i) = x0(i) - SUM( invers_H(i,1:n)*g(1:n) )
          end do
          line_search = .false.
       else
          if (outp > 1) write (*,*) 'Newton_Opt_2D: problem concave, line search',f
          ft_min = f
          length_scale = 0.1 * dot_product( x0, x0 )
          direction(1:n) = g(1:n) / dot_product( g, g )
          xt(1:n) = x0(1:n) - length_scale * direction(1:n)
          do i = 1, 2*line_steps
             xt(1:n) = xt(1:n) + length_scale * direction(1:n) / real( line_steps )
             call f_val ( ft, xt, n )
             if ( ft < ft_min ) then
                x(1:n) = xt(1:n)
                ft_min = ft
             end if
          end do
          if (outp > 1) write (*,*) 'finished line search', x(1:n)
          line_search = .true.
       end if

       call f_val (ftrial, x, n)

       if ( ftrial > f .AND. .NOT. line_search ) then
          if (outp > 1) write (*,*) 'reduced step size'
          do i = 1,n
             x(i) = x0(i) - 0.5 * SUM( invers_H(i,1:n)*g(1:n) )
          end do
          ! call f_val (ftrial, x, n)
          ! if ( ftrial > f ) then
          !   do i = 1,n
          !      x(i) = x0(i) - 0.25 * SUM( invers_H(i,1:n)*g(1:n) )
          !   end do
          ! end if
       end if

       error_in_x = SUM( ABS( x(1:n)-x0(1:n) ) )

       if (outp > 1) write (*,'(a,4(G15.7))') ' finished Newton step: f,err,par ',f,SUM( ABS( g(1:n) ) ), x(1:n)

       ! write (*,*) 'f= ',fmin
       ! write (*,*) 'par',exp( x(1) ), exp( x(2) )
       ! write (*,'(a,4(G15.7))') ' I= ',invers_H(1,1),invers_H(1,2),invers_H(2,1),invers_H(2,2)
       ! write (*,'(a,4(G15.7))') ' H= ',hessian(1,1),hessian(1,2),hessian(2,1),hessian(2,2)

    end do


  end subroutine Newton_Opt_2D



  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  SUBROUTINE f_grad_hessian (hessian, f_val, gtrans, fmin, optpara, n)

    !USE BASIC_VARIABLES
    !USE STARTING_VALUES
    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: n
    REAL, INTENT(IN)           :: optpara(:)
    REAL, INTENT(IN OUT)       :: fmin
    REAL, INTENT(IN OUT)       :: gtrans(:)
    REAL, INTENT(IN OUT)       :: hessian(:,:)

    INTERFACE
       SUBROUTINE f_val (f, x, n)
         REAL, INTENT(IN OUT)                :: f
         REAL, INTENT(IN)                    :: x(:)
         INTEGER, INTENT(IN)                 :: n
       END SUBROUTINE f_val
       !SUBROUTINE f_grad (g, x, n)
       !  REAL, INTENT(IN OUT)                     :: g(:)
       !  REAL, INTENT(IN)                         :: x(:)
       !  INTEGER, INTENT(IN)                      :: n
       !END SUBROUTINE f_grad
    END INTERFACE

    !---------------------------------------------------------------------------
    INTEGER                                 :: i, j
    REAL                                    :: delta
    REAL                                    :: optpara_mod(n), g(n), gi(n,n), gi_left(n,n)
    !---------------------------------------------------------------------------


    delta = 2.0E-5

    optpara_mod = optpara

    DO i = 1, n

       optpara_mod = optpara
       optpara_mod(i) = optpara(i)*(1.0+delta)

       CALL f_grad (g, f_val, optpara_mod, n)
       gi(i,1:n) = g(1:n)

       optpara_mod(i) = optpara(i)*(1.0-delta)

       CALL f_grad (g, f_val, optpara_mod, n)
       gi_left(i,1:n) = g(1:n)

    ENDDO

    CALL f_grad (g, f_val, optpara, n)


    DO i = 1, n
       DO j = 1, n

          hessian(i,j) = ( gi(i,j) - gi_left(i,j) ) / ( 2.0*optpara(i)*delta )
          ! hessian(j,i) = hessian(i,j)
          ! write (*,*) i,j,hessian(i,j)

       ENDDO
    ENDDO

    gtrans = g
    CALL f_val (fmin, optpara, n)

  END SUBROUTINE f_grad_hessian


  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  !
  !WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  SUBROUTINE f_grad (g, f_val, optpara, n)

    !USE basic_variables
    !USE starting_values
    IMPLICIT NONE

    !-----------------------------------------------------------------------------
    INTEGER, INTENT(IN)                    :: n
    REAL, INTENT(IN)                       :: optpara(:)
    REAL, INTENT(IN OUT)                   :: g(:)

    INTERFACE
       SUBROUTINE f_val ( fmin, optpars, n )
         INTEGER, INTENT(IN)                :: n
         REAL, INTENT(IN)                   :: optpars(:)
         REAL, INTENT(IN OUT)               :: fmin
       END SUBROUTINE f_val
    END INTERFACE

    !-----------------------------------------------------------------------------
    INTEGER                                :: i
    REAL                                   :: delta, fmin
    REAL                                   :: optpara_mod(n), fi(n), f
    !-----------------------------------------------------------------------------


    delta = 1.E-5

    optpara_mod = optpara

    DO i = 1, n

       optpara_mod = optpara
       optpara_mod(i) = optpara(i)*(1.0+delta)

       call f_val ( fmin, optpara_mod, n )
       fi(i) = fmin

    ENDDO

    call f_val ( fmin, optpara, n )
    f = fmin


    DO i = 1, n

       g(i) = ( fi(i) - f ) / ( optpara(i)*delta )
       !  write (*,*) fi(i), f, optpara(i)*delta

    ENDDO

  END SUBROUTINE f_grad


END MODULE optimizer_2D
