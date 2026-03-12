!> @mainpage Main program
!> @authors
!> Rishabh Bhatt, Laurent Debreu, Arthur Vidard

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Program: main
!
! Description:
!   Solves the linear system obtained when using the Parareal algorithm for the 
!   tangent-linear integrations in the inner loop of incremental 4D-Var. The 
!   model here is the linearised 2D rotating shallow water equations and the 
!   minimisation uses the Conjugate Gradient (CG) method
!   Three CG variants are tested:
!     1. Exact CG with exact matrix
!     2. Exact CG with Parareal matrix approximation
!     3. Inexact CG with Parareal
!
! Dependencies:
!   - shallow_water
!   - parareal_utils
!   - eigenvalues
!   - output_netcdf
!   - shallow_water_adjoint
!   - inexact_version
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

program main

use shallow_water
use parareal_utils
use eigenvalues
use output_netcdf
use shallow_water_adjoint
use inexact_version

integer :: i
integer :: it, nitermax

real :: t1,t2,t3,t4,t5,t6,t7,t8
real :: rnorm,xout,eps1

real, dimension(:), allocatable :: obs,B,x0,x1,x2,test,r,y
real, dimension(:,:), allocatable :: cg_solution, obs_sp  ! obs_sp: array of sparse observations

double precision :: dnrm2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initialise shallow water model
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

call initialise_sw()
call CreatenetCDF()

! Compute the true solution using sequential fine solver
allocate(Truesolution(FineGrid%nx_1D,0:N_time_windows))
Truesolution(:,0) = initial_xn

call OutPutNetCDF()

write(*,*)
!call my_cpu_time(t3)
do i=1,N_time_windows

  write(*,*) 'Time window i = ',i, 'max(eta)', maxval(FineGrid%eta)
  call integrate_sw(etan=FineGrid%eta,un=Finegrid%u,vn=FineGrid%v,nt=Nfine)

  call trans_2Dto1D(FineGrid%u,FineGrid%v,FineGrid%eta,Truesolution(:,i),FineGrid%nx,FineGrid%ny, &
                FineGrid%nx_1D,FineGrid%indices_1D)
  call OutPutNetCDF()
end do
call my_cpu_time(t4)
!t4=t4-t3
!write(*,*),'Timing non parareal = ',t4
write(*,*)

call CloseNetCDF()


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Parareal initialisation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write(*, *) repeat('#', 30)

call initialise_parareal(FineGrid,CoarseGrid)
write(*,*) 'Parareal initialised'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Observation initialisation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allocate(B,x0,x1,x2,y,mold=initial_xn)
allocate(obs_sp(FineGrid%nx_1D,Nobs))

call initialise_obs_mask(FineGrid%nx_1D)
call initialise_obs(initial_xn,FineGrid%nx_1D,obs_sp)

call Bvector(FineGrid%nx_1D,obs_sp,B)
write(*,*) 'B vector computed'

x0 = 0.

write(*,*) 'Begin conjgrad'
call set_verbose(.TRUE.)

if (regularisation) then
  write(*,*) 'Using regularisation constant, alpharegul = ', alpharegul
else
  write(*,*) 'No regularisation used'
end if

write (*,*) repeat('#',10), 'Test 1: Run exact CG with exact matrix ', repeat('#',10)
write (*,*)
!call my_cpu_time(t1)
call conjgrad(Mmatrix,B,x0,FineGrid%nx_1D,eps=1.,observation=obs_sp,exact_matrix=Mmatrix)
!call my_cpu_time(t2)
!write(*,*) 'timing CG = ',t2-t1 


x1 = 0.
write (*,*)
write(*, *) repeat('#', 10), ' Test 2: Exact CG with Parareal matrix ', repeat('#', 10)
write (*,*)
!call my_cpu_time(t5)
call conjgrad(Mmatrix_parareal,B,x1,FineGrid%nx_1D,eps=1.D0,observation=obs_sp,exact_matrix=Mmatrix)
!call my_cpu_time(t6)
!write(*,*) 'timing CG = ',t6-t5
!stop


! Compute tolerance for inexact CG
allocate(r(FineGrid%nx_1D),test(FineGrid%nx_1D))

call Mmatrix(FineGrid%nx_1D,x0,test)
r = test - B

call matrixnorminv(Mmatrix,r,FineGrid%nx_1D,rnorm,eps=1.)
write(*,*) 'Residual norm (rnorm) = ', rnorm

call matrixnorminv(Mmatrix,B,FineGrid%nx_1D,xout,eps=1.)
write(*,*) 'RHS norm (bnorm) = ', xout

call get_tol(rnorm,xout,eps1)

write(*,*) "inexact CG tolerance (eps) = ",eps1
deallocate(r)

x2 = 0.
nitermax=20
write (*,*)
write(*, *) repeat('#', 10), ' Test 3: Inexact CG with Parareal ', repeat('#', 10)
!call my_cpu_time(t7)
call inexact_conjgrad(B,x2,FineGrid%nx_1D,eps=eps1,nitermax=nitermax,observation=obs_sp)
!call my_cpu_time(t8)
!write(*,*) 'timing inexact CG = ',t8-t7

call give_output(cg_solution)


contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Subroutine: give_output
  !
  ! Description:
  !   Runs the forward shallow water model from the CG-retrieved initial
  !   condition (x1) and writes the resulting trajectory to a NetCDF file
  !   named 'cg_output.nc'.
  !
  ! Arguments:
  !   sol (out) - 2D array storing the model state at each time window
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine give_output(sol)
  
  real, dimension(:,:), allocatable :: sol
  
  ! to change the name of file to cg_output.nc, else overwrites output.nc
  call change_output(.true.)
  
  call initialise_sw(x1)    ! uses x1 to initialise the shallow water model, initial state can be changed
  call CreateNetCDF()
  allocate(sol(FineGrid%nx_1D,0:N_time_windows))
  
  initial_xn = x1
  sol(:,0) = initial_xn 
  
  call OutPutNetCDF()
  
  do i=1,N_time_windows
    write(*,*) 'Time window i = ',i, 'max(eta)', maxval(FineGrid%eta)
    call integrate_sw(etan=FineGrid%eta,un=FineGrid%u,vn=FineGrid%v,nt=Nfine)
  
    call trans_2Dto1D(FineGrid%u,FineGrid%v,FineGrid%eta,sol(:,i),FineGrid%nx,FineGrid%ny, &
                  FineGrid%nx_1D,FineGrid%indices_1D)
    call outPutNetCDF()
  end do
  
  call CloseNetCDF()
  deallocate(sol)
  
  end subroutine give_output  




end program main
