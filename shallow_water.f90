!> Shallow water module

!> @brief
!! Utilities to run the shallow water model

module shallow_water
use shallow_water_utils
use shallow_water_step

real :: Lx      !< Length of the domain
real :: Ly      !< Length of the domain
real :: H       !< Depth of the domain

real :: g       !< gravity

real :: coriolis !< Coriolis Parameter

real,private :: solver_prec   !< solver precision
integer, parameter :: max_iterations = 1000 !< maximum iterations for gmres (used for inverting the matrix)


!> @brief create discretisation grids
!> @details gives grids of different step size which in turn are used to construct the fine and coarse solvers
type grid
real,dimension(:,:),allocatable :: u   !< velocity x-component
real,dimension(:,:),allocatable :: v   !< velocity y-component
real,dimension(:,:),allocatable :: eta !< free surface elevation
real,dimension(:,:),allocatable :: xr  !< cells center x-location
real,dimension(:,:),allocatable :: yr  !< cells center y-location
real,dimension(:),allocatable :: max_eigenvalue_sw  !< maximum eigenvalue (in amplitude) of the SW matrix
real :: dt     !< time step
real :: dx     !< grid step x-diection
real :: dy     !< grid step y-direction
integer :: nx  !< Number of cells in the domain, x-direction
integer :: ny  !< Number of cells in the domain, y-direction
integer :: nx_1D  !< size of the vector containing velocities and free surface
integer,dimension(:,:,:),allocatable :: indices_1D  !< array of indices, the third dimension is for variables

integer :: nt_time_windows !< Number of time steps per time windows
end type grid

type(grid),pointer :: FineGrid  !< Fine grid
type(grid),pointer :: CoarseGrid  !< Coarse grid

type(grid),pointer :: curgrid !< Pointer on the currently integrated grid

real :: omegamax !< Courant number
real :: theta    !< coefficient of the theta scheme

real,dimension(:,:),allocatable :: initial_eta  !< initial free surface
real,dimension(:,:),allocatable :: initial_u    !< initial velocity in x-direction
real,dimension(:,:),allocatable :: initial_v    !< initial velocity in y-direction  
real,dimension(:),allocatable :: initial_xn     !< initial 1D array containing initial-eta, initial_u, and initial_v

real,dimension(:,:),allocatable :: Truesolution !< True solution (for testing purpose)

contains

real function get_solver_prec()
  get_solver_prec = solver_prec
end function get_solver_prec

subroutine set_solver_prec(my_prec)
  real :: my_prec
  solver_prec = my_prec
end subroutine set_solver_prec


subroutine initialise_sw(init1)
 use eigenvalues

real, dimension(curgrid%nx_1D), optional :: init1

 solver_prec=1.D-12
 allocate(FineGrid,CoarseGrid)

! Initialize Coarse grid : dx, nx_1D, indices_1D
 curgrid => CoarseGrid

 call initialise_sw_parameters()
 call compute_nx_1D(curgrid%nx,curgrid%ny,curgrid%nx_1D,curgrid%indices_1D)

! Initialize Fine grid : dx, nx_1D, indices_1D
 curgrid => FineGrid

 call initialise_sw_parameters()
 call compute_nx_1D(curgrid%nx,curgrid%ny,curgrid%nx_1D,curgrid%indices_1D)

if (present(init1)) then
 call initialise_sw_variables(init1)
else
 call initialise_sw_variables()
end if

! compute the largest eigenvalue of the RHS (matrix associated to the SW equations)
 call compute_eigenvalues(shallow_water_matrix,curgrid%nx_1D, &
                          curgrid%max_eigenvalue_sw,svd=.false.)

 print *,'max_eigenvalue_sw = ',curgrid%max_eigenvalue_sw

 curgrid%dt = omegamax / curgrid%max_eigenvalue_sw(1)

 !curgrid%dt=curgrid%dt
 print *,'Fine time step dt = ',curgrid%dt

end subroutine initialise_sw

!> Initialize the shallow water parameters
subroutine initialise_sw_parameters

Lx = 200000.
Ly = 200000.
H = 100.
!g = 9.81
!g=0.066
g=2.

!coriolis = 1.e-2
coriolis = 5.e-3

curgrid%nx = 80
curgrid%ny = 80
curgrid%dx = Lx/curgrid%nx
curgrid%dy = Lx/curgrid%ny

omegamax = 1.1  !< Courant number, for finding a reasonable time-step
theta = 0.51    !< theta parameter for the implicit numerical finite difference scheme

print *,'Rossby Radius = ',sqrt(g*H)/coriolis

end subroutine initialise_sw_parameters

subroutine set_theta(theta_in)
  real :: theta_in
  
  theta=theta_in
end subroutine set_theta


!> @brief Initialise shallow water variables
!> @details variables - horizontal velocities u and v, free surface eta
!! - allocate appropriate storage
!! - define grid point values
!! - put initial values to u,v and eta using variables initial_u, initial_v, and initial_eta
!! - transforming the 2D array of initial values (eta,u,v) to a 1D array initial_xn
subroutine initialise_sw_variables(init)
integer :: i,j
integer :: nx, ny
real, dimension(curgrid%nx_1D), optional :: init

  real,dimension(1:curgrid%nx+1,0:curgrid%ny+1) :: un
  real,dimension(0:curgrid%nx+1,1:curgrid%ny+1) :: vn
  real,dimension(0:curgrid%nx+1,0:curgrid%ny+1) :: etan
  real, dimension(1:curgrid%nx+1,0:curgrid%ny+1) :: rhs_u
  real, dimension(0:curgrid%nx+1,1:curgrid%ny+1) :: rhs_v
  real, dimension(0:curgrid%nx+1,0:curgrid%ny+1) :: rhs_eta

nx = curgrid%nx
ny = curgrid%ny

if (allocated(curgrid%u)) then
  deallocate(curgrid%u)
end if

if (allocated(curgrid%v)) then
  deallocate(curgrid%v)
end if

if (allocated(curgrid%eta)) then
  deallocate(curgrid%eta)
end if

if (allocated(curgrid%xr)) then
  deallocate(curgrid%xr)
end if

if (allocated(curgrid%yr)) then
  deallocate(curgrid%yr)
end if


allocate(curgrid%u(1:nx+1,0:ny+1))
allocate(curgrid%v(0:nx+1,1:ny+1))
allocate(curgrid%eta(0:nx+1,0:ny+1))
allocate(curgrid%xr(0:nx+1,0:ny+1))
allocate(curgrid%yr(0:nx+1,0:ny+1))


do j=0,ny+1
do i=0,nx+1
  curgrid%xr(i,j) = i*curgrid%dx-curgrid%dx/2.
  curgrid%yr(i,j) = j*curgrid%dy-curgrid%dy/2.
enddo
enddo

curgrid%u = 0.
curgrid%v = 0.

if (present(init)) then
  call trans_1Dto2D(init,un,vn,etan,nx,ny,curgrid%nx_1D,curgrid%indices_1D)
  curgrid%eta = etan
  curgrid%u = un
  curgrid%v = vn
else
! assign a gaussian function at the middle of the basin
do j=0,ny+1
do i=0,nx+1
  curgrid%eta(i,j) = exp(-((curgrid%xr(i,j)-0.5*Lx)/(Lx/10.))**2-((curgrid%yr(i,j)-Ly/2.)/(Ly/10.))**2)
enddo
enddo
end if

if (allocated(initial_u)) then
  deallocate(initial_u)
end if

if (allocated(initial_v)) then
  deallocate(initial_v)
end if

if (allocated(initial_eta)) then
  deallocate(initial_eta)
end if

if (allocated(initial_xn)) then
  deallocate(initial_xn)
end if


allocate(initial_u(1:nx+1,0:ny+1))
allocate(initial_v(0:nx+1,1:ny+1))
allocate(initial_eta(0:nx+1,0:ny+1))
allocate(initial_xn(curgrid%nx_1D))

initial_u = curgrid%u
initial_v = curgrid%v
initial_eta=curgrid%eta

call trans_2Dto1D(initial_u,initial_v,initial_eta,initial_xn,nx,ny,curgrid%nx_1D,curgrid%indices_1D)

allocate(curgrid%max_eigenvalue_sw(1))

end subroutine initialise_sw_variables

!> Integrate the shallow water model for a number of time step nt
!> @param[in] nt
subroutine integrate_sw(xn,etan,un,vn,nt)
  use eigenvalues
  real,dimension(curgrid%nx_1D),optional :: xn
  real,dimension(0:curgrid%nx+1,0:curgrid%ny+1),optional :: etan
  real,dimension(1:curgrid%nx+1,0:curgrid%ny+1),optional :: un
  real,dimension(0:curgrid%nx+1,1:curgrid%ny+1),optional :: vn

  integer :: nt !> number of time steps

  integer :: nb,i
  real,dimension(curgrid%nx_1D) :: rhs_xn,xcur
  real,parameter :: eps = 1.D-8

  integer :: nx,ny,nx_1D
  real :: dt

!$ integer :: omp_get_thread_num

  nx = curgrid%nx
  ny = curgrid%ny
  dt = curgrid%dt
  nx_1D = curgrid%nx_1D

!  print *,'dt = ',curgrid%dt,nt
  if (present(etan)) then
    call trans_2Dto1D(un,vn,etan,xcur,nx,ny,curgrid%nx_1D,curgrid%indices_1D)
  else
    xcur = xn
  endif

  do nb=1,nt
    !print *,'nb = ',nb
    call shallow_water_matrix(nx_1D,xcur,rhs_xn)

      rhs_xn = xcur + dt*(1.-theta)*rhs_xn
   call gmres ( shallow_water_theta_matrix, nx_1D, xcur, rhs_xn, &
          itr_max=max_iterations, mr=max(10,nx_1D/100), tol_abs=solver_prec, tol_rel=solver_prec)
  enddo
  
  if (present(etan)) then
    call trans_1Dto2D(xcur,un,vn,etan,nx,ny,curgrid%nx_1D,curgrid%indices_1D)
  else
    xn = xcur
  endif

end subroutine integrate_sw

subroutine compute_matrix_sw()
  real,dimension(curgrid%nx_1D) :: xn,yn


  integer :: i,j
  real,dimension(curgrid%nx_1D) :: rhs_xn
  real :: val
  integer :: nx_1D


  nx_1D = curgrid%nx_1D

  do j=1,nx_1D
    yn=0.
    yn(j)=1.
  do i=1,nx_1D
    xn=0.
    xn(i)=1.
    call shallow_water_matrix(nx_1D,xn,rhs_xn)
    val = DOT_PRODUCT(yn,rhs_xn)
  enddo
enddo
end subroutine compute_matrix_sw


!> @brief gives the matrix-vector product of the right hand side of the shallow water matrix
subroutine shallow_water_matrix(n,xn,rhs_xn,k)
  integer :: n
  real,dimension(n) :: xn
  real,dimension(n) :: rhs_xn
  integer, optional :: k

  real,dimension(1:curgrid%nx+1,0:curgrid%ny+1) :: un
  real,dimension(0:curgrid%nx+1,1:curgrid%ny+1) :: vn
  real,dimension(0:curgrid%nx+1,0:curgrid%ny+1) :: etan
  real, dimension(1:curgrid%nx+1,0:curgrid%ny+1) :: rhs_u
  real, dimension(0:curgrid%nx+1,1:curgrid%ny+1) :: rhs_v
  real, dimension(0:curgrid%nx+1,0:curgrid%ny+1) :: rhs_eta


   ! 1D array xn of length curgrid%nx_1D contains values of eta,u, and v
   ! returns 2D arrays of etan,un, and vn 
    call trans_1Dto2D(xn,un,vn,etan,curgrid%nx,curgrid%ny,curgrid%nx_1D,curgrid%indices_1D)

   ! uses values of 2D arrays etan,un, and vn for computing the discretised rhs of shallow water
   ! returns 2D arrays rhs_eta,rhs_u, and rhs_v 
    call rhs_sw(un,vn,etan,rhs_u,rhs_v,rhs_eta,curgrid%nx,curgrid%ny,g,H,coriolis,curgrid%dx,curgrid%dy)

   ! returns a 1D array rhs_xn containing rhs_eta,rhs_u, and rhs_v
   ! rhs_xn --> gives A*x
    call trans_2Dto1D(rhs_u,rhs_v,rhs_eta,rhs_xn,curgrid%nx,curgrid%ny,curgrid%nx_1D,curgrid%indices_1D)

end subroutine shallow_water_matrix


!> @brief gives the matrix-vector product of the shallow water matrix when theta scheme is used

!> If we write the shallow water equations as a system i.e. \f$ \mathrm{d} X/ \mathrm{d} t = AX \f$ \n
!> using implicit theta scheme we have \f$ (X^{n+1} - X^n)/ \delta t = \theta A X^{n+1} + (1-\theta) A X^n \f$ \n
!> simplifying leads to \f$[I- \delta t \theta A] X^{n+1} = [I+ \delta t (1-\theta) A] X^n \f$ \n
!> \n
!> The subroutine returns \f$ [I- \delta t \theta A] X \f$
subroutine shallow_water_theta_matrix(n,xn,rhs_xn,k)
  integer :: n
  real,dimension(n) :: xn
  real,dimension(n) :: rhs_xn
  integer, optional :: k

  call shallow_water_matrix(n,xn,rhs_xn)

  rhs_xn = xn-curgrid%dt*theta*rhs_xn

end subroutine shallow_water_theta_matrix

end module shallow_water
