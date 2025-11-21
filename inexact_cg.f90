!> implements the inexact conjugate gradient
module inexact_version

use parareal_utils

contains



!> @brief compute the true residual i.e. \f$ \underbrace{(F^N)^T(F^N)}_{A}\textbf{x} - \underbrace{(F^N)^T\textbf{y}}_{\textbf{b}} \f$

!> @param[in] nx size
!> @param x vector
!> @param bvec the right hand side vector
subroutine residual(nx,x,bvec,y)

  integer :: nx
  real, dimension(nx) :: x, bvec
  real, dimension(nx) :: y

  real, dimension(nx) :: y1,y2
  
  y=x
  call Mmatrix(nx,y,y1)

  y = y1- bvec

end subroutine residual


subroutine cost_func(nx,x,y,ob)

  integer :: nx,i
  real, dimension(nx) :: x,y1,ob
  real :: y
  double precision :: dnrm2
  
  real :: res

  call Fmatrix(nx,x,y1)
  y = 0.5*dnrm2(nx,y1-ob,1)

  if (regularisation) then
    y = y + 0.5*alpharegul*(dnrm2(nx,x,1))**2
  end if

end subroutine cost_func


subroutine get_tol(rnrm,bnrm,e1) 

  real :: rnrm, bnrm
  real, intent(out) :: e1

  e1 = (2*rnrm/bnrm)**2

end subroutine get_tol



subroutine gradient(nx,x,res,bvec)
  
  integer :: nx
  double precision :: dnrm2
  
  real, dimension(nx) :: res1, x, bvec, y1

  real :: res

  call Mmatrix(nx,x,res1)

  res = dnrm2(nx,res1 - bvec,1)
  
  if (regularisation) then
    call identity_matrix(nx,x,y1)
    res = res + dnrm2(nx,alpharegul*y1,1)
  end if

end subroutine gradient


subroutine quadratic(nx,x,quad,bvec)
  
  integer :: nx
  real, dimension(nx) :: x, Ax, bvec

  real :: quad

  call Mmatrix(nx,x,Ax)
  quad = 0.5*DOT_PRODUCT(x,Ax) - DOT_PRODUCT(bvec,x)

end subroutine quadratic 


!> @brief computes the ellipsoidal norm of a matrix.
!> \f$ \Vert x \Vert_A = x^TAx \f$
subroutine ep_norm(matrix_v,x,n,res)
  
  procedure(matrix_vector) :: matrix_v
  integer :: n
  real, dimension(n) :: x, Ax

  real, intent(out) :: res

  call matrix_v(n,x,Ax)

  res = SQRT(DOT_PRODUCT(x,Ax))

end subroutine ep_norm


!> @brief gives the stopping tolerance for the primal-dual norm of the error matrix
!> \f$ \Vert E_j \Vert_{A^{-1},A} \leq \omega_j \f$

!> @param eps \f$ \epsilon \f$
!> @param bnorm \f$ \Vert b \Vert_{A^{-1}} \f$
!> @param pnorm \f$ \Vert p \Vert_A \f$
!> @param rnorm \f$ \Vert r \Vert_2 \f$
!> @param phi \f$ \phi \f$
subroutine E_bound(eps,bnorm,pnorm,rnorm,phi,xout)
  
  real :: eps, phi, bnorm, pnorm, rnorm
  real, intent(out) :: xout

  xout = (SQRT(eps)*bnorm*pnorm)/(2*phi*(rnorm**2) + SQRT(eps)*bnorm*pnorm)
  
end subroutine E_bound


subroutine inexact_quadratic(nx,x,bvec,q)

  integer :: nx
  real, dimension(nx) :: x, bvec
  real :: q

  q = -0.5*DOT_PRODUCT(bvec,x)

end subroutine inexact_quadratic



subroutine inexact_conjgrad(b,x,n,eps,nitermax,observation)
  
  integer :: n
  real, dimension(n) :: b, x, x_new, Ax, c, c1, c3
  real, dimension(n) :: r, p, p_new, r_new
  real, dimension(:), allocatable :: q1, cf1
  
  real :: eps, phi, phi_new, omega, omega2, e, e1, c2, c4
  real :: alpha, beta_old, beta_new, bound
  real :: omega_hat, phi_hat, rinvnorm

  integer :: iteration, nitermax
  real :: big_phi_new, big_phi = 1.
  integer :: d, k, i, l
  real, dimension(:,:), allocatable :: para_sol, p_sol, p_sol1 

  logical :: reorth = .true.
  logical :: inacc_budget = .true.
  logical :: cp

  integer :: k1, it_out, ksol, k_bar, k2
  real :: xout, xout1, xout2, q, cf, q2
  real :: bnorm, pnorm, r2norm, exact_pnorm, exact_bnorm, exact_ebound
  real, dimension(n) :: tmp, tmp1, tmp3, tmp4, tmp5, x_est, obs, gap1, gap2 
  real, dimension(:), allocatable :: res, sol1, sol2
  real, dimension(:,:), allocatable :: u
  real :: res1, res2, tmp2, xi1, xi2
  
  real :: eps1

  double precision :: dnrm2
  logical :: exists

  real, dimension(0:nitermax) :: quad_vals
  
  real, dimension(n), optional :: observation
 
  if (regularisation) then
    print *, 'regularisation constant for inexact cg', alpharegul
  end if

  print *, 'size ', n 
  beta_old = (dnrm2(n,b,1))**2
  print *,'beta old = ', beta_old 

  !call matrixnorminv(Mmatrix,b,n,exact_bnorm)
  !exact_bnorm = 9.6549753709801731  

  !alpharegul = 5
  exact_bnorm = 9.6557121601659652 
  print *, 'bnorm = ', exact_bnorm

  eps1 = 1.      

  r = -b
  p = b

  bound = 0.5*SQRT(eps)*exact_bnorm
  print *, 'Stopping tolerance = ',bound

  iteration = 0
  phi = nitermax
  print *, 'nitermax = ', nitermax

  allocate(u(n,0:2*N_time_windows))
  u(:,0) = b/beta_old
  
  quad_vals(0) = 0

  d = 2 

  do while (.true.)

    print *, '***********************************************'
    print *, 'Inexact CG Iteration: ', iteration
    
    if (iteration .gt. 0) then 

      inacc_budget = .true.
    
      print *, 'phi = ', phi
    
      bnorm = SQRT(2*ABS(q2))
      !bnorm = exact_bnorm
      print *, 'bnorm = ', bnorm

      call ep_norm(Mmatrix,p,n,xout)
      exact_pnorm = xout

      print *, 'exact pnorm = ', exact_pnorm

      r2norm = dnrm2(n,r,1)
      print *, 'r2norm = ', r2norm

      call E_bound(eps,exact_bnorm,exact_pnorm,r2norm,phi,xout)
      omega = xout

      print *, 'omega = ', omega
    
      omega2 = omega*exact_pnorm
      print *, 'omega2 = ', omega2

    end if 

    inquire(file="exact_norm.csv",exist = exists)
    if(exists) then
       open(4,file='exact_norm.csv',form='formatted',status = 'old', position = 'append',action='write')
    else
       open(4,file='exact_norm.csv',form='formatted',status = 'new', position = 'append',action='write')
    end if

    write(4,'(F15.10,F15.10,F15.10,F15.10)') exact_pnorm, omega, omega2, bnorm  
    close(4)

    !!!!!! with  exact ||Ep||_A^{-1} !!!!!!
    !k1 = 1 
    
    !if(allocated(para_sol)) then
      !deallocate(para_sol)
    !end if

    !call set_pverbose(.false.)
    !call parareal(FineGrid,CoarseGrid,p,sol1,maxk=8,iter_out=it_out,full_sol=para_sol)

    !call FMatrix(n,p,tmp1) 

    !do while (k1<=N_time_windows)
      
      !e = dnrm2(n,para_sol(:,k1) - tmp1,1)
      !print *, 'exact ep norm at parareal iteration ', k1, ' = ', e

      !if (e > omega2) then
        !k1 = k1+1
        !cycle
      !else
        !omega_hat = e
        !print *, 'omega hat: ', omega_hat
        !ksol = k1
        !print *, 'Number of parareal iterations = ', ksol

        !c = para_sol(:,ksol)
        !do i =N_time_windows-1,0,-1
          !call INTEGRATE_SW_ADJOINT(xnb=c,nt=curgrid%nt_time_windows)
        !end do
        !if (regularisation) then
          !call identity_matrix(n,p,tmp3)
          !c = c + alpharegul*tmp3
        !end if
        !exit
      !end if
    !end do


!!!!!!!!!!!! Using last parareal iterate !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.true.) then
    if (iteration .eq. 0) then
      call Mmatrix(n,p,c)
    end if


    if (iteration .gt. 0) then

    if (allocated(para_sol)) then
      deallocate(para_sol)
    end if
    
    call set_pverbose(.false.)

    call parareal(FineGrid,CoarseGrid,p,sol1,maxk=10,iter_out=it_out,full_sol=para_sol)
    print *, 'no of iter check', it_out
    
    write(*,*)

    do k1 = 1,N_time_windows

      !!!!!!!!! p norm approximations !!!!!!!!!
      c1 = para_sol(:,k1)
      c2 = DOT_PRODUCT(c1,c1)

      write(*,*)
      if (regularisation) then
        call identity_matrix(n,p,tmp4)
        pnorm = sqrt(c2 + alpharegul*DOT_PRODUCT(p,tmp4))
      else
        pnorm = sqrt(c2)
      end if
      print *, 'approx pnorm at iteration ', k1, '=', pnorm

      call E_bound(eps,exact_bnorm,pnorm,r2norm,phi,xi1)
      xout1 = xi1*pnorm 
      print *, 'omega2 with approx pnorm ',xout1

      call E_bound(eps,bnorm,exact_pnorm,r2norm,phi,xi2)
      xout2 = xi2*exact_pnorm
      print *, 'omega2 with approx bnorm ',xout2

      call E_bound(eps,bnorm,pnorm,r2norm,phi,xout)
      omega = xout

      print *, 'omega = ', omega
    
      omega2 = omega*pnorm
      print *, 'omega2 = ', omega2

      e = dnrm2(n,para_sol(:,k1+1)-para_sol(:,k1),1)
      print *, 'approx ep norm at iteration ', k1, ' = ', e

    inquire(file="para_norms.csv",exist = exists)
    if(exists) then
       open(5,file='para_norms.csv',form='formatted',status = 'old', position = 'append',action='write')
    else
       open(5,file='para_norms.csv',form='formatted',status = 'new', position = 'append',action='write')
    end if

    write(5,'(F15.10,F15.10,F15.10,F15.10,F15.10,F15.10)') pnorm, omega, omega2, e, xout1, xout2
    close(5)

      write(*,*)
      if (e < omega2) then
        omega_hat = e
        print *, 'omega_hat : ', omega_hat
        print *, 'parareal iterations ', k1+1

        c3 = para_sol(:,k1+1)
        c4 = DOT_PRODUCT(c3,c3)

        if (regularisation) then
          call identity_matrix(n,p,tmp5)
          pnorm = sqrt(c4 + alpharegul*DOT_PRODUCT(p,tmp5))
        else
          pnorm = sqrt(c4)
        end if
        print *, 'actual pnorm used at iteration ', k1+1, '=', pnorm
        ksol = k1
        exit
      end if
    end do

    write(*,*)

    call FMatrix(n,p,tmp1) 
    do k2 =1,ksol
      e1 = dnrm2(n,para_sol(:,k2) - tmp1,1)
      print *, 'exact ep norm at parareal iteration ', k2, ' = ', e1

      inquire(file="exact_ep.csv",exist = exists)
      if(exists) then
        open(10,file='exact_ep.csv',form='formatted',status = 'old', position = 'append',action='write')
      else
        open(10,file='exact_ep.csv',form='formatted',status = 'new', position = 'append',action='write')
      end if

      write(10,'(F15.10)') e1  
      close(10)

    end do
    
    print *, 'ksol for c calculation', ksol

    c = para_sol(:,ksol)
    do i =N_time_windows-1,0,-1
      call INTEGRATE_SW_ADJOINT(xnb=c,nt=curgrid%nt_time_windows)
    end do
    if (regularisation) then
      call identity_matrix(n,p,tmp3)
      c = c + alpharegul*tmp3
    end if

    if (inacc_budget) then
      phi_hat = ((pnorm - omega_hat)*SQRT(eps)*bnorm*pnorm)/(omega_hat*2*((r2norm)**2))
      !phi_hat = ((exact_pnorm - omega_hat)*SQRT(eps)*exact_bnorm*exact_pnorm)/(omega_hat*2*((r2norm)**2))

      print *, 'phi: ', phi, 'phi_hat: ', phi_hat, 'phi_hat > phi :', phi_hat>phi

      big_phi_new = big_phi - (1/phi_hat)

      if (iteration < nitermax) then
        phi_new = (nitermax-iteration-1)/big_phi_new
      else
        phi_new = phi
      end if

      phi = phi_new
      big_phi = big_phi_new

    end if
    end if
    
    end if


    alpha = beta_old/DOT_PRODUCT(p,c)
    x_new = x + alpha*p
    r_new = r + alpha*c
    
    if (iteration .gt. 0) then
      
      write(*,*)
      if(present(observation)) then
        allocate(cf1(n))
        call Fmatrix(n,x_new,cf1)
        cf = 0.5*(dnrm2(n,cf1-observation,1))**2
        if (regularisation) then
          cf = cf + 0.5*alpharegul*(dnrm2(n,x,1))**2
        end if
        print *, 'cost func value: ', cf
        deallocate(cf1)
      end if

      call inexact_quadratic(n,x_new,b,q2)
      quad_vals(iteration+1) = q2
      print *, 'approx quadratic at iteration ', iteration, ' = ', q2

      allocate(q1(n))
      call Mmatrix(n,x_new,q1)
      q = 0.5*DOT_PRODUCT(x_new,q1) - DOT_PRODUCT(b,x_new)
      print *, 'exact quadratic value ', q
      deallocate(q1)
      
      print *, '2-norm of r :', dnrm2(n,r_new,1)
    
      call matrixnorminv(Mmatrix,r_new,n,xout,eps=eps1)
      rinvnorm = xout

      print *, 'A-1 norm of r :',rinvnorm 

      call Mmatrix(n,x_new,gap1)
      gap2 = gap1 - b
      call matrixnorminv(Mmatrix,gap2 - r_new,n,res2,eps=eps1)
      print *, 'residual gap norm: ', res2 
      write(*,*)

    else
      call inexact_quadratic(n,x_new,b,q2)
      print *, 'approx quadratic at iteration ', iteration, ' = ', q2
      quad_vals(iteration+1) = q2 

      allocate(q1(n))
      call Mmatrix(n,x_new,q1)
      q = 0.5*DOT_PRODUCT(x_new,q1) - DOT_PRODUCT(b,x_new)
      print *, 'exact quadratic value ', q
      deallocate(q1)
      
    end if

    inquire(file="icg_norm.csv",exist = exists)
    if(exists) then
       open(3,file='icg_norm.csv',form='formatted',status = 'old', position = 'append',action='write')
    else
       open(3,file='icg_norm.csv',form='formatted',status = 'new', position = 'append',action='write')
    end if

    write(3,'(F15.10,F15.10,F15.10,F15.10,F15.10,F15.10,F15.10)') sqrt(DOT_PRODUCT(r_new,r_new)),rinvnorm,omega,omega2,cf,q,res2
    close(3)

    ! stopping criterion
    !if (rinvnorm <= bound) exit

    !!!! approximate stopping criterion !!!!
    if (iteration >=d) then
      print *, 'quad diff :', quad_vals(iteration+1-d) - quad_vals(iteration+1)
      if ( (quad_vals(iteration+1-d) - quad_vals(iteration+1)) .le. 0.25*eps*ABS(quad_vals(iteration+1)) ) exit
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (reorth) then

      do i = 0,iteration
        r_new = r_new - DOT_PRODUCT(u(:,i),r_new)*u(:,i)
      end do

      beta_new = DOT_PRODUCT(r_new,r_new)
      u(:,iteration+1) = r_new/SQRT(beta_new)
    else
      beta_new = DOT_PRODUCT(r_new,r_new)

    end if

    p_new = - r_new + (beta_new/beta_old)*p

    ! Update old values

    beta_old = beta_new
    r = r_new
    x = x_new
    p = p_new

    iteration = iteration +1
  end do

print *, 'CG iterations: ', iteration

end subroutine inexact_conjgrad


end module inexact_version
