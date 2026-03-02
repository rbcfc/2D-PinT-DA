!> implements the inexact conjugate gradient
module inexact_version

use parareal_utils
use eigenvalues

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
  real, dimension(nx,Nobs) :: y1, ob
  real, dimension(nx) :: x
  real :: y, temp
  double precision :: dnrm2
  
  real :: res

  call Fmatrix(nx,x,y1)

  y=0.
  do i = 1,Nobs
    temp = 0.5*dnrm2(nx,y1(:,i)-ob(:,i),1)
    y = y + temp
  end do

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
  real, dimension(n) :: b, x, x_new, c, c1, ep_k, ep_k1
  real, dimension(n) :: r, p, p_new, r_new
  real, dimension(:), allocatable :: q1
  real, dimension(:,:), allocatable :: cf1
  
  real :: eps, phi, phi_new, omega, omega2, e, e1, c2, c4, c5, c6
  real :: alpha, beta_old, beta_new, bound
  real :: omega_hat, phi_hat, rinvnorm

  integer :: iteration, nitermax
  real :: big_phi_new, big_phi = 1.
  integer :: d, k, i, j 
  real, dimension(:,:), allocatable :: para_sol

  logical :: reorth = .true.
  logical :: inacc_budget = .true.

  integer :: k1, it_out, ksol, k2, icf
  real :: xout, xout1, xout2, q, cf, q2, cf_temp
  real :: bnorm, pnorm, r2norm, exact_pnorm, exact_bnorm, exact_ebound
  real, dimension(n) :: temp, tmp1, tmp3, tmp4, tmp5, tmp6, gap1, gap2 
  real, dimension(n) :: ep_temp, ep_temp1
  real, dimension(:), allocatable :: sol1
  real, dimension(:,:), allocatable :: u
  real :: res2, xi1, xi2
  
  integer, dimension(nitermax) :: it_array
  real, dimension(n,N_time_windows,N_time_windows) :: all_iter

  double precision :: dnrm2
  logical :: file_exists

  real, dimension(0:nitermax) :: quad_vals
  
  real, dimension(n,Nobs), optional :: observation
 
  if (regularisation) then
    print *, 'regularisation constant for inexact cg', alpharegul
  end if

  print *, 'size ', n 
  beta_old = (dnrm2(n,b,1))**2
  print *,'beta old = ', beta_old 

  !if (Nobs > 1) then
  !  call matrixnorminv(Mmatrix,b,n,exact_bnorm,eps=1.)
  !else
  !  call matrixnorminv(Mmatrix,b,n,exact_bnorm,eps=1.D-1)
  !end if

  !exact_bnorm = 9.6557121601659652 
  exact_bnorm = 11.894898781675595

  print *, 'bnorm = ', exact_bnorm

  r = -b
  p = b

  bound = 0.5*SQRT(eps)*exact_bnorm
  print *, 'Stopping tolerance = ',bound

  iteration = 0
  phi = nitermax
  print *, 'nitermax = ', nitermax

  allocate(q1(n),u(n,0:2*N_time_windows))
  u(:,0) = b/beta_old
  
  quad_vals(0) = 0

  d = 1 

  inquire(file="inexact_cg.csv", exist=file_exists)

  if (.not. file_exists) then
    open(3, file='inexact_cg.csv', form='formatted', status='new')
    write(3, '(A)') 'iteration,r2norm,rinvnorm,omega,omega2,cost_func,quadratic,res_gap'
    close(3)
  end if

  cf = 0.0
  q = 0.0
  res2 = 0.0
  rinvnorm = 0.0
  omega = 0.0
  omega2 = 0.0

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

    inquire(file="exact_norm.csv",exist = file_exists)
    if(file_exists) then
       open(4,file='exact_norm.csv',form='formatted',status = 'old', position = 'append',action='write')
    else
       open(4,file='exact_norm.csv',form='formatted',status = 'new', position = 'append',action='write')
    end if

    write(4,'(F15.10,F15.10,F15.10,F15.10)') exact_pnorm, omega, omega2, bnorm  
    close(4)


    !!!!!!!!!!!! Using last parareal iterate for ||Ep||_A^-1 approximation !!!!!!!!!!!!
    if (iteration .eq. 0) then
      call Mmatrix(n,p,c)
    end if

    if (iteration .gt. 0) then

      if (allocated(para_sol)) then
        deallocate(para_sol)
      end if
    
      call set_pverbose(.false.)

      call parareal(FineGrid,CoarseGrid,p,sol1,maxk=13,iter_out=it_out,full_sol=para_sol,all_iterates=all_iter)
      print *, 'no of iter check', it_out

      write(*,*)
      
      do k1 = 1,N_time_windows

        !!!!!!!!! p norm approximations !!!!!!!!!
        c2 = 0.0
        do i = 1, Nobs
          c1 = all_iter(:,obs_windows(i),k1)
          where (.not. obs_mask(:, i)) c1 = 0.0

          c6 = DOT_PRODUCT(c1,c1)
          c2 = c2+c6
        end do

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

        
        ! ep norm approximation
        ep_k = 0.
        ep_k1 = 0.

        do i = 1, Nobs
          ep_temp = all_iter(:,obs_windows(i),k1)
          ep_temp1 = all_iter(:,obs_windows(i),k1+1)
          where (.not. obs_mask(:, i)) ep_temp  = 0.0
          where (.not. obs_mask(:, i)) ep_temp1 = 0.0

          ep_k = ep_k+ep_temp
          ep_k1 = ep_k1+ep_temp1
        end do

        e = dnrm2(n,ep_k-ep_k1,1)
        print *, 'approx ep norm at iteration ', k1, ' = ', e

        inquire(file="para_norms.csv",exist = file_exists)
        if(file_exists) then
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
          
          c4 = 0.0
          do i = 1, Nobs
            temp = all_iter(:,obs_windows(i),k1+1)
            where (.not. obs_mask(:, i)) temp = 0.0
            c5 = DOT_PRODUCT(temp,temp)
            c4 = c4 +c5
          end do

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

      call Fmatrix(n,p,tmp1) 
      do k2 =1,ksol
        e1 = dnrm2(n,para_sol(:,k2) - tmp1,1)
        print *, 'exact ep norm at parareal iteration ', k2, ' = ', e1

        inquire(file="exact_ep.csv",exist = file_exists)
        if(file_exists) then
          open(10,file='exact_ep.csv',form='formatted',status = 'old', position = 'append',action='write')
        else
          open(10,file='exact_ep.csv',form='formatted',status = 'new', position = 'append',action='write')
        end if

        write(10,'(F15.10)') e1  
        close(10)

      end do
    
      print *, 'ksol for c calculation', ksol

      c = 0.
      do i = 1, Nobs
        tmp6 = all_iter(:,obs_windows(i),ksol)
        where (.not. obs_mask(:, i)) tmp6 = 0.0

        do j = obs_windows(i)-1,0,-1
          call INTEGRATE_SW_ADJOINT(xnb=tmp6,nt=curgrid%nt_time_windows)
        end do
        c = c + tmp6
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


    alpha = beta_old/DOT_PRODUCT(p,c)
    x_new = x + alpha*p
    r_new = r + alpha*c
    
    if (iteration .gt. 0) then
      
      write(*,*)
      if(present(observation)) then
        allocate(cf1(n,Nobs))
        call Fmatrix(n,x_new,cf1)
        cf = 0.0
        do icf = 1,Nobs
          cf_temp = 0.5*(dnrm2(n,cf1(:,icf) - observation(:,icf),1))**2
          cf = cf + cf_temp
        end do
        if (regularisation) then
          cf = cf + 0.5*alpharegul*(dnrm2(n,x_new,1))**2
        end if
        print *, 'cost func value: ', cf
        deallocate(cf1)
      end if

      call inexact_quadratic(n,x_new,b,q2)
      quad_vals(iteration+1) = q2
      print *, 'approx quadratic at iteration ', iteration, ' = ', q2

      call Mmatrix(n,x_new,q1)
      q = 0.5*DOT_PRODUCT(x_new,q1) - DOT_PRODUCT(b,x_new)
      print *, 'exact quadratic value ', q
      
      print *, '2-norm of r :', dnrm2(n,r_new,1)
    
      call matrixnorminv(Mmatrix,r_new,n,xout,eps=1.)
      rinvnorm = xout

      !print *, 'A-1 norm of r :',rinvnorm 

      call Mmatrix(n,x_new,gap1)
      gap2 = gap1 - b
      call matrixnorminv(Mmatrix,gap2 - r_new,n,res2,eps=1.)
      print *, 'residual gap norm: ', res2 
      write(*,*)

    else
      call inexact_quadratic(n,x_new,b,q2)
      print *, 'approx quadratic at iteration ', iteration, ' = ', q2
      quad_vals(iteration+1) = q2 

      call Mmatrix(n,x_new,q1)
      q = 0.5*DOT_PRODUCT(x_new,q1) - DOT_PRODUCT(b,x_new)
      print *, 'exact quadratic value ', q
      
    end if

    open(3, file='inexact_cg.csv', form='formatted', status='old', position='append')
    write(3, '(I5,A,F15.10,A,F15.10,A,F15.10,A,F15.10,A,F15.10,A,F15.10,A,F15.10)') &
      iteration, ',', sqrt(DOT_PRODUCT(r_new,r_new)), ',', rinvnorm, ',', &
      omega, ',', omega2, ',', cf, ',', q, ',', res2
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

  deallocate(q1,u)

print *, 'CG iterations: ', iteration

end subroutine inexact_conjgrad


subroutine conjgrad(matrix_v,b,x,n,eps,observation,exact_matrix)


  procedure(matrix_vector) :: matrix_v
  procedure(matrix_vector), optional :: exact_matrix
  integer :: n,i,icf
  real,dimension(n) :: b,x,obs
  real :: eps

  real,dimension(n) :: r,Ax,Ap,p,Apold,pold,ptmp
  integer :: iteration

  real :: rsold,rsnew,alpha,tmp,cf,res,tmp2,grad,quad
  real, dimension(:), allocatable :: Aptmp
  real, dimension(:,:), allocatable :: cftmp

  logical :: file_exists
  double precision :: dnrm2

  real, dimension(n,Nobs), optional :: observation
  real :: equad, cf_temp

  call matrix_v(n,x,Ax)

  r = b -Ax

  if (verbose_cg .eqv. .true.) then
    print *,'r = ',DOT_PRODUCT(r,r)
    print *,'bnorm = ',DOT_PRODUCT(b,b),sqrt(DOT_PRODUCT(b,b))
  end if
  p = r
  rsold = DOT_PRODUCT(r,r)

  equad = -46.491068930860784

  iteration = 1

  inquire(file="conjgrad.csv", exist=file_exists)

  if (.not. file_exists) then
    open(2, file='conjgrad.csv', form='formatted', status='new')
    write(2, '(A)') 'iteration,r2norm,cost_func,quadratic,Ainv_rnorm'
    close(2)
  end if

  cf   = 0.0
  quad = 0.0
  tmp  = 0.0

  do while (.TRUE.)

    print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print *, 'CG iteration number ', iteration
    write (*,*)
    call matrix_v(n,p,Ap)

    alpha = rsold / DOT_PRODUCT(p,Ap)
    x = x + alpha * p
    r = r - alpha * Ap
    rsnew = DOT_PRODUCT(r,r)

    if (verbose_cg .eqv. .true.) then
      print *,'ratio = ',rsnew/rsold
      print *,'iter = ',iteration,rsnew,'r2-norm: ',sqrt(rsnew)
    end if
   
     if(present(observation)) then
       allocate(cftmp(n,Nobs))
       call Fmatrix(n,x,cftmp)
       cf = 0.0
       do icf = 1,Nobs
         cf_temp = 0.5*(dnrm2(n,cftmp(:,icf) - observation(:,icf),1))**2
         cf = cf + cf_temp
       end do
       if (regularisation) then
         cf = cf + 0.5*alpharegul*(dnrm2(n,x,1))**2
       end if
       print *, 'cost func value: ', cf
       deallocate(cftmp)
     end if
    !call matrixnorminv(exact_matrix,r,n,tmp)
    !print *, 'A-1 norm of r: ',tmp
    
    if (PRESENT(exact_matrix)) then
      allocate(Aptmp(n))
      call exact_matrix(n,x,Aptmp)

      quad = 0.5*DOT_PRODUCT(x,APtmp) - DOT_PRODUCT(b,x)
      print *, 'quadratic value : ', quad

      tmp =  SQRT(2*(quad - equad))
      print *, 'A-1 norm of r', tmp
      deallocate(Aptmp)
    end if

    open(2, file='conjgrad.csv', form='formatted', status='old', position='append')
    write(2, '(I5,A,F15.10,A,F15.10,A,F15.10,A,F15.10)') &
      iteration, ',', sqrt(rsnew), ',', cf, ',', quad, ',', tmp
    close(2)

    if (sqrt(rsnew) < eps) exit

    p = r + (rsnew/rsold) * p
    rsold = rsnew

    iteration = iteration + 1
    !print *, 'interation counter changes to iter = ', iteration
  enddo

  if (verbose_cg .eqv. .true.) then
    print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print *,'lastiter = ',iteration,rsnew,sqrt(rsnew)
  end if

end subroutine conjgrad

end module inexact_version
